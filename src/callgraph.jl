@enum Direction FFT_FORWARD=-1 FFT_BACKWARD=1
@enum Pow24 POW2=2 POW4=1
@enum FFTEnum compositeFFT dft pow2FFT pow3FFT pow4FFT bluesteinFFT

"""
$(TYPEDSIGNATURES)
Node of a call graph

# Arguments
`left`: Offset to the left child node
`right`: Offset to the right child node
`type`: Object representing the type of FFT
`sz`: Size of this FFT

"""
struct CallGraphNode{T}
    left::Int
    right::Int
    type::FFTEnum
    sz::Int
    s_in::Int
    s_out::Int
    w::T
end

"""
$(TYPEDSIGNATURES)
Object representing a graph of FFT Calls

# Arguments
`nodes`: Nodes keeping track of the graph
`workspace`: Preallocated Workspace

"""
struct CallGraph{T<:Complex}
    nodes::Vector{CallGraphNode{T}}
    workspace::Vector{Vector{T}}
end

# Get the node in the graph at index i
Base.getindex(g::CallGraph{T}, i::Int) where {T} = g.nodes[i]

"""
$(TYPEDSIGNATURES)
Check if `N` is a power of 2 or 4

"""
function _ispow24(N::Int)
    N < 1 && return nothing
    while N & 0b11 == 0
        N >>= 2
    end
    return N < 3 ? Pow24(N) : nothing
end

"""
$(TYPEDSIGNATURES)
Recursively instantiate a set of `CallGraphNode`s

# Arguments
`nodes`: A vector (which gets expanded) of `CallGraphNode`s
`N`: The size of the FFT
`workspace`: A vector (which gets expanded) of preallocated workspaces
`s_in`: The stride of the input
`s_out`: The stride of the output
`d`: The direction of the transform

"""
function CallGraphNode!(nodes::Vector{CallGraphNode{T}}, N::Int, workspace::Vector{Vector{T}}, s_in::Int, s_out::Int, d::Direction)::Int where {T}
    if N == 0
        throw(DimensionMismatch("array has to be non-empty"))
    end
    w = cispi(T(2)/N)
    if iseven(N)
        pow = _ispow24(N)
        if !isnothing(pow)
            push!(workspace, T[])
            push!(nodes, CallGraphNode(0, 0, pow == POW2 ? pow2FFT : pow4FFT, N, s_in, s_out, w))
            return 1
        end
    end
    if N % 3 == 0
        if nextpow(3, N) == N
            push!(workspace, T[])
            push!(nodes, CallGraphNode(0, 0, pow3FFT, N, s_in, s_out, w))
            return 1
        end
    end
    if N == 1 || Primes.isprime(N)
        if N >= 13
            # Allocate and precompute workspace for Bluestein algorithm
            # Workspace layout: [chirp(N), b_fft(M), a(M), a_fft(M)]
            # where M = nextpow(2, 2*N-1)
            # chirp is precomputed here, b_fft will be lazily computed on first use
            M = nextpow(2, 2*N - 1)
            tmp = zeros(T, 2*M + N)

            # Precompute chirp sequence: chirp[n] = w^(n²/2)
            # Apply direction: for forward use conj(w), for backward use w
            w_dir = d == FFT_FORWARD ? conj(w) : w
            w_half = sqrt(w_dir)
            chirp_power = one(T)
            chirp_mult = w_half
            @inbounds for n in 0:N-1
                tmp[n+1] = chirp_power
                chirp_power *= chirp_mult
                chirp_mult *= w_dir
            end

            # Store b vector (not yet FFT'd) in the b_fft section for now
            # It will be transformed to FFT(b) on first call to fft_bluestein!
            @inbounds for n in 0:N-1
                tmp[N+n+1] = conj(tmp[n+1])  # b[n] = conj(chirp[n])
            end
            @inbounds for n in 1:N-1
                tmp[N+M-n+1] = conj(tmp[n+1])  # b[M-n] = conj(chirp[n])
            end

            push!(workspace, tmp)
            push!(nodes, CallGraphNode(0, 0, bluesteinFFT, N, s_in, s_out, w))
        else
            push!(workspace, T[])
            push!(nodes, CallGraphNode(0, 0, dft, N, s_in, s_out, w))
        end
        return 1
    end
    Ns = [first(x) for x in collect(Primes.factor(N)) for _ in 1:last(x)]
    if Ns[1] == 2
        N1 = prod(Ns[Ns .== 2])
    elseif Ns[1] == 3
        N1 = prod(Ns[Ns .== 3])
    else
        # Greedy search for closest factor of N to sqrt(N)
        Nsqrt = sqrt(N)
        N_cp = cumprod(Ns[end:-1:1])[end:-1:1]
        N_prox = abs.(N_cp .- Nsqrt)
        _,N1_idx = findmin(N_prox)
        N1 = N_cp[N1_idx]
    end
    N2 = N ÷ N1
    push!(nodes, CallGraphNode(0, 0, dft, N, s_in, s_out, w))
    sz = length(nodes)
    push!(workspace, Vector{T}(undef, N))
    left_len = CallGraphNode!(nodes, N1, workspace, N2, N2*s_out, d)
    right_len = CallGraphNode!(nodes, N2, workspace, N1*s_in, 1, d)
    nodes[sz] = CallGraphNode(1, 1 + left_len, compositeFFT, N, s_in, s_out, w)
    return 1 + left_len + right_len
end

"""
$(TYPEDSIGNATURES)
Instantiate a CallGraph from a number `N`

"""
function CallGraph{T}(N::Int, d::Direction=FFT_FORWARD) where {T}
    nodes = CallGraphNode{T}[]
    workspace = Vector{Vector{T}}()
    CallGraphNode!(nodes, N, workspace, 1, 1, d)
    CallGraph(nodes, workspace)
end
