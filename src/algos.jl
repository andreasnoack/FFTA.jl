@inline function direction_sign(d::Direction)
    Int(d)
end

@inline _conj(w::Complex, d::Direction) = ifelse(direction_sign(d) === 1, w, conj(w))

function fft!(out::AbstractVector{T}, in::AbstractVector{T}, start_out::Int, start_in::Int, d::Direction, t::FFTEnum, g::CallGraph{T}, idx::Int) where T
    if t === compositeFFT
        fft_composite!(out, in, start_out, start_in, d, g, idx)
    else
        root = g[idx]
        if t == dft
            fft_dft!(out, in, root.sz, start_out, root.s_out, start_in, root.s_in, _conj(root.w, d))
        else
            N = root.sz
            s_in = root.s_in
            s_out = root.s_out
            if t === pow2FFT
                fft_pow2!(out, in, N, start_out, s_out, start_in, s_in, _conj(root.w, d))
            elseif t === pow3FFT
                p_120 = cispi(T(2)/3)
                m_120 = cispi(T(4)/3)
                _p_120, _m_120 = d == FFT_FORWARD ? (p_120, m_120) : (m_120, p_120)
                fft_pow3!(out, in, N, start_out, s_out, start_in, s_in, _conj(root.w, d), _m_120, _p_120)
            elseif t === pow4FFT
                fft_pow4!(out, in, N, start_out, s_out, start_in, s_in, _conj(root.w, d))
            elseif t === pow8FFT
                fft_pow8!(out, in, N, start_out, s_out, start_in, s_in, _conj(root.w, d))
            else
                throw(ArgumentError("kernel not implemented"))
            end
        end
    end
end


"""
$(TYPEDSIGNATURES)
Cooley-Tukey composite FFT, with a pre-computed call graph

# Arguments
`out`: Output vector
`in`: Input vector
`start_out`: Index of the first element of the output vector
`start_in`: Index of the first element of the input vector
`d`: Direction of the transform
`g`: Call graph for this transform
`idx`: Index of the current transform in the call graph

"""
function fft_composite!(out::AbstractVector{T}, in::AbstractVector{U}, start_out::Int, start_in::Int, d::Direction, g::CallGraph{T}, idx::Int) where {T,U}
    root = g[idx]
    left_idx = idx + root.left
    right_idx = idx + root.right
    left = g[left_idx]
    right = g[right_idx]
    N  = root.sz
    N1 = left.sz
    N2 = right.sz
    s_in = root.s_in
    s_out = root.s_out

    w1 = _conj(root.w, d)
    wj1 = one(T)
    tmp = g.workspace[idx]
    @inbounds for j1 in 0:N1-1
        wk2 = wj1
        fft!(tmp, in, N2*j1+1, start_in + j1*s_in, d, right.type, g, right_idx)
        j1 > 0 && @inbounds for k2 in 1:N2-1
            tmp[N2*j1 + k2 + 1] *= wk2
            wk2 *= wj1
        end
        wj1 *= w1
    end

    @inbounds for k2 in 0:N2-1
        fft!(out, tmp, start_out + k2*s_out, k2+1, d, left.type, g, left_idx)
    end
end

"""
$(TYPEDSIGNATURES)
Discrete Fourier Transform, O(N^2) algorithm, in place.

# Arguments
`out`: Output vector
`in`: Input vector
`N`: Size of the transform
`start_out`: Index of the first element of the output vector
`stride_out`: Stride of the output vector
`start_in`: Index of the first element of the input vector
`stride_in`: Stride of the input vector
`w`: The value `cispi(direction_sign(d) * 2 / N)`

"""
function fft_dft!(out::AbstractVector{T}, in::AbstractVector{T}, N::Int, start_out::Int, stride_out::Int, start_in::Int, stride_in::Int, w::T) where {T}
    tmp = in[start_in]
    @inbounds for j in 1:N-1
        tmp += in[start_in + j*stride_in]
    end
    out[start_out] = tmp

    wk = wkn = w
    @inbounds for d in 1:N-1
        tmp = in[start_in]
        @inbounds for k in 1:N-1
            tmp += wkn*in[start_in + k*stride_in]
            wkn *= wk
        end
        out[start_out + d*stride_out] = tmp
        wk *= w
        wkn = wk
    end
end

function fft_dft!(out::AbstractVector{Complex{T}}, in::AbstractVector{T}, N::Int, start_out::Int, stride_out::Int, start_in::Int, stride_in::Int, w::Complex{T}) where {T<:Real}
    halfN = N÷2

    tmp = Complex{T}(in[start_in])
    @inbounds for j in 1:N-1
        tmp += in[start_in + j*stride_in]
    end
    out[start_out] = tmp

    wk = wkn = w
    @inbounds for d in 1:halfN
        tmp = Complex{T}(in[start_in])
        @inbounds for k in 1:N-1
            tmp += wkn*in[start_in + k*stride_in]
            wkn *= wk
        end
        out[start_out + d*stride_out] = tmp
        wk *= w
        wkn = wk
    end
end

"""
$(TYPEDSIGNATURES)
Power of 2 FFT, in place

# Arguments
`out`: Output vector
`in`: Input vector
`N`: Size of the transform
`start_out`: Index of the first element of the output vector
`stride_out`: Stride of the output vector
`start_in`: Index of the first element of the input vector
`stride_in`: Stride of the input vector
`w`: The value `cispi(direction_sign(d) * 2 / N)`

"""
function fft_pow2!(out::AbstractVector{T}, in::AbstractVector{U}, N::Int, start_out::Int, stride_out::Int, start_in::Int, stride_in::Int, w::T) where {T, U}
    if N == 2
        out[start_out]              = in[start_in] + in[start_in + stride_in]
        out[start_out + stride_out] = in[start_in] - in[start_in + stride_in]
        return
    end
    m = N ÷ 2

    fft_pow2!(out, in, m, start_out               , stride_out, start_in            , stride_in*2, w*w)
    fft_pow2!(out, in, m, start_out + m*stride_out, stride_out, start_in + stride_in, stride_in*2, w*w)

    wj = one(T)
    @inbounds for j in 0:m-1
        j1_out = start_out + j*stride_out
        j2_out = start_out + (j+m)*stride_out
        out_j    = out[j1_out]
        out[j1_out] = out_j + wj*out[j2_out]
        out[j2_out] = out_j - wj*out[j2_out]
        wj *= w
    end
end


"""
$(TYPEDSIGNATURES)
Power of 4 FFT, in place

# Arguments
`out`: Output vector
`in`: Input vector
`N`: Size of the transform
`start_out`: Index of the first element of the output vector
`stride_out`: Stride of the output vector
`start_in`: Index of the first element of the input vector
`stride_in`: Stride of the input vector
`w`: The value `cispi(direction_sign(d) * 2 / N)`

"""
function fft_pow4!(out::AbstractVector{T}, in::AbstractVector{U}, N::Int, start_out::Int, stride_out::Int, start_in::Int, stride_in::Int, w::T) where {T, U}
    minusi = -sign(imag(w))*im
    @inbounds if N == 4
        xee = in[start_in]
        xoe = in[start_in +   stride_in]
        xeo = in[start_in + 2*stride_in]
        xoo = in[start_in + 3*stride_in]
        xee_p_xeo = xee + xeo
        xee_m_xeo = xee - xeo
        xoe_p_xoo = xoe + xoo
        xoe_m_xoo = -(xoe - xoo)*minusi
        out[start_out]                = xee_p_xeo + xoe_p_xoo
        out[start_out +   stride_out] = xee_m_xeo + xoe_m_xoo
        out[start_out + 2*stride_out] = xee_p_xeo - xoe_p_xoo
        out[start_out + 3*stride_out] = xee_m_xeo - xoe_m_xoo
        return
    end
    m = N ÷ 4

    w1 = w
    w2 = w*w1
    w3 = w*w2
    w4 = w*w3

    fft_pow4!(out, in, m, start_out                 , stride_out, start_in              , stride_in*4, w4)
    fft_pow4!(out, in, m, start_out +   m*stride_out, stride_out, start_in +   stride_in, stride_in*4, w4)
    fft_pow4!(out, in, m, start_out + 2*m*stride_out, stride_out, start_in + 2*stride_in, stride_in*4, w4)
    fft_pow4!(out, in, m, start_out + 3*m*stride_out, stride_out, start_in + 3*stride_in, stride_in*4, w4)

    wkoe = wkeo = wkoo = one(T)

    @inbounds for k in 0:m-1
        kee = start_out +  k          * stride_out
        koe = start_out + (k +     m) * stride_out
        keo = start_out + (k + 2 * m) * stride_out
        koo = start_out + (k + 3 * m) * stride_out
        y_kee, y_koe, y_keo, y_koo = out[kee], out[koe], out[keo], out[koo]
        ỹ_keo = y_keo*wkeo
        ỹ_koe = y_koe*wkoe
        ỹ_koo = y_koo*wkoo
        y_kee_p_y_keo = y_kee + ỹ_keo
        y_kee_m_y_keo = y_kee - ỹ_keo
        ỹ_koe_p_ỹ_koo = ỹ_koe + ỹ_koo
        ỹ_koe_m_ỹ_koo = -(ỹ_koe - ỹ_koo) * minusi
        out[kee] = y_kee_p_y_keo + ỹ_koe_p_ỹ_koo
        out[koe] = y_kee_m_y_keo + ỹ_koe_m_ỹ_koo
        out[keo] = y_kee_p_y_keo - ỹ_koe_p_ỹ_koo
        out[koo] = y_kee_m_y_keo - ỹ_koe_m_ỹ_koo
        wkoe *= w1
        wkeo *= w2
        wkoo *= w3
    end
end


"""
$(TYPEDSIGNATURES)
Power of 8 FFT, in place

# Arguments
`out`: Output vector
`in`: Input vector
`N`: Size of the transform
`start_out`: Index of the first element of the output vector
`stride_out`: Stride of the output vector
`start_in`: Index of the first element of the input vector
`stride_in`: Stride of the input vector
`w`: The value `cispi(direction_sign(d) * 2 / N)`

"""
function fft_pow8!(out::AbstractVector{T}, in::AbstractVector{U}, N::Int, start_out::Int, stride_out::Int, start_in::Int, stride_in::Int, w::T) where {T, U}
    minusi = -sign(imag(w))*im
    @inbounds if N == 8
        # Base case: 8-point FFT decomposed as two 4-point FFTs
        # Even indices: 0, 2, 4, 6
        # Odd indices: 1, 3, 5, 7
        x0 = in[start_in]
        x1 = in[start_in +   stride_in]
        x2 = in[start_in + 2*stride_in]
        x3 = in[start_in + 3*stride_in]
        x4 = in[start_in + 4*stride_in]
        x5 = in[start_in + 5*stride_in]
        x6 = in[start_in + 6*stride_in]
        x7 = in[start_in + 7*stride_in]

        # 4-point FFT of even indices (0,2,4,6)
        e0e2 = x0 + x4
        e0_2 = x0 - x4
        e2e6 = x2 + x6
        e2_6 = (x2 - x6) * minusi
        E0 = e0e2 + e2e6
        E1 = e0_2 + e2_6
        E2 = e0e2 - e2e6
        E3 = e0_2 - e2_6

        # 4-point FFT of odd indices (1,3,5,7)
        o1o5 = x1 + x5
        o1_5 = x1 - x5
        o3o7 = x3 + x7
        o3_7 = (x3 - x7) * minusi
        O0 = o1o5 + o3o7
        O1 = o1_5 + o3_7
        O2 = o1o5 - o3o7
        O3 = o1_5 - o3_7

        # Twiddle factors for combining
        # W = e^(-2πi/8) = e^(-πi/4)
        w8_1 = cispi(T(-1)/4)  # e^(-πi/4)
        w8_2 = minusi          # e^(-πi/2) = -i (for forward)
        w8_3 = w8_1 * w8_2     # e^(-3πi/4)

        # Combine: X[k] = E[k] + W^k * O[k], X[k+4] = E[k] - W^k * O[k]
        out[start_out]                = E0 + O0
        out[start_out +   stride_out] = E1 + w8_1 * O1
        out[start_out + 2*stride_out] = E2 + w8_2 * O2
        out[start_out + 3*stride_out] = E3 + w8_3 * O3
        out[start_out + 4*stride_out] = E0 - O0
        out[start_out + 5*stride_out] = E1 - w8_1 * O1
        out[start_out + 6*stride_out] = E2 - w8_2 * O2
        out[start_out + 7*stride_out] = E3 - w8_3 * O3
        return
    end

    m = N ÷ 8

    w1 = w
    w2 = w*w1
    w3 = w*w2
    w4 = w*w3
    w5 = w*w4
    w6 = w*w5
    w7 = w*w6
    w8 = w*w7

    # Recursively process 8 subproblems
    fft_pow8!(out, in, m, start_out                 , stride_out, start_in              , stride_in*8, w8)
    fft_pow8!(out, in, m, start_out +   m*stride_out, stride_out, start_in +   stride_in, stride_in*8, w8)
    fft_pow8!(out, in, m, start_out + 2*m*stride_out, stride_out, start_in + 2*stride_in, stride_in*8, w8)
    fft_pow8!(out, in, m, start_out + 3*m*stride_out, stride_out, start_in + 3*stride_in, stride_in*8, w8)
    fft_pow8!(out, in, m, start_out + 4*m*stride_out, stride_out, start_in + 4*stride_in, stride_in*8, w8)
    fft_pow8!(out, in, m, start_out + 5*m*stride_out, stride_out, start_in + 5*stride_in, stride_in*8, w8)
    fft_pow8!(out, in, m, start_out + 6*m*stride_out, stride_out, start_in + 6*stride_in, stride_in*8, w8)
    fft_pow8!(out, in, m, start_out + 7*m*stride_out, stride_out, start_in + 7*stride_in, stride_in*8, w8)

    # Combine with 8-point butterflies following the same decomposition pattern
    wk1 = wk2 = wk3 = wk4 = wk5 = wk6 = wk7 = one(T)
    w8_twiddle1 = cispi(T(-1)/4)  # W_8^1
    w8_twiddle2 = minusi          # W_8^2
    w8_twiddle3 = w8_twiddle1 * w8_twiddle2  # W_8^3

    @inbounds for k in 0:m-1
        k0 = start_out +  k          * stride_out
        k1 = start_out + (k +     m) * stride_out
        k2 = start_out + (k + 2 * m) * stride_out
        k3 = start_out + (k + 3 * m) * stride_out
        k4 = start_out + (k + 4 * m) * stride_out
        k5 = start_out + (k + 5 * m) * stride_out
        k6 = start_out + (k + 6 * m) * stride_out
        k7 = start_out + (k + 7 * m) * stride_out

        y0, y1, y2, y3 = out[k0], out[k1], out[k2], out[k3]
        y4, y5, y6, y7 = out[k4], out[k5], out[k6], out[k7]

        # Apply twiddle factors from the recursive calls
        ỹ1 = y1 * wk1
        ỹ2 = y2 * wk2
        ỹ3 = y3 * wk3
        ỹ4 = y4 * wk4
        ỹ5 = y5 * wk5
        ỹ6 = y6 * wk6
        ỹ7 = y7 * wk7

        # 4-point FFT of even-indexed results (0,2,4,6)
        e0e4 = y0 + ỹ4
        e0_4 = y0 - ỹ4
        e2e6 = ỹ2 + ỹ6
        e2_6 = (ỹ2 - ỹ6) * minusi
        E0 = e0e4 + e2e6
        E1 = e0_4 + e2_6
        E2 = e0e4 - e2e6
        E3 = e0_4 - e2_6

        # 4-point FFT of odd-indexed results (1,3,5,7)
        o1o5 = ỹ1 + ỹ5
        o1_5 = ỹ1 - ỹ5
        o3o7 = ỹ3 + ỹ7
        o3_7 = (ỹ3 - ỹ7) * minusi
        O0 = o1o5 + o3o7
        O1 = o1_5 + o3_7
        O2 = o1o5 - o3o7
        O3 = o1_5 - o3_7

        # Combine: X[k] = E[k] + W_8^k * O[k], X[k+4] = E[k] - W_8^k * O[k]
        out[k0] = E0 + O0
        out[k1] = E1 + w8_twiddle1 * O1
        out[k2] = E2 + w8_twiddle2 * O2
        out[k3] = E3 + w8_twiddle3 * O3
        out[k4] = E0 - O0
        out[k5] = E1 - w8_twiddle1 * O1
        out[k6] = E2 - w8_twiddle2 * O2
        out[k7] = E3 - w8_twiddle3 * O3

        wk1 *= w1
        wk2 *= w2
        wk3 *= w3
        wk4 *= w4
        wk5 *= w5
        wk6 *= w6
        wk7 *= w7
    end
end


"""
$(TYPEDSIGNATURES)
Power of 3 FFT, in place

# Arguments
out: Output vector
in: Input vector
N: Size of the transform
start_out: Index of the first element of the output vector
stride_out: Stride of the output vector
start_in: Index of the first element of the input vector
stride_in: Stride of the input vector
w: The value `cispi(direction_sign(d) * 2 / N)`
plus120: Depending on direction, perform either ±120° rotation
minus120: Depending on direction, perform either ∓120° rotation

"""
function fft_pow3!(out::AbstractVector{T}, in::AbstractVector{U}, N::Int, start_out::Int, stride_out::Int, start_in::Int, stride_in::Int, w::T, plus120::T, minus120::T) where {T, U}
    if N == 3
        @muladd out[start_out + 0]            = in[start_in] + in[start_in + stride_in]          + in[start_in + 2*stride_in]
        @muladd out[start_out +   stride_out] = in[start_in] + in[start_in + stride_in]*plus120  + in[start_in + 2*stride_in]*minus120
        @muladd out[start_out + 2*stride_out] = in[start_in] + in[start_in + stride_in]*minus120 + in[start_in + 2*stride_in]*plus120
        return
    end

    # Size of subproblem
    Nprime = N ÷ 3

    # Dividing into subproblems
    fft_pow3!(out, in, Nprime, start_out, stride_out, start_in, stride_in*3, w^3, plus120, minus120)
    fft_pow3!(out, in, Nprime, start_out + Nprime*stride_out, stride_out, start_in + stride_in, stride_in*3, w^3, plus120, minus120)
    fft_pow3!(out, in, Nprime, start_out + 2*Nprime*stride_out, stride_out, start_in + 2*stride_in, stride_in*3, w^3, plus120, minus120)

    w1 = w
    w2 = w*w1
    wk1 = wk2 = one(T)
    for k in 0:Nprime-1
        @muladd k0 = start_out + k*stride_out
        @muladd k1 = start_out + (k+Nprime)*stride_out
        @muladd k2 = start_out + (k+2*Nprime)*stride_out
        y_k0, y_k1, y_k2 = out[k0], out[k1], out[k2]
        @muladd out[k0] = y_k0 + y_k1*wk1 + y_k2*wk2
        @muladd out[k1] = y_k0 + y_k1*wk1*plus120 + y_k2*wk2*minus120
        @muladd out[k2] = y_k0 + y_k1*wk1*minus120 + y_k2*wk2*plus120
        wk1 *= w1
        wk2 *= w2
    end
end
