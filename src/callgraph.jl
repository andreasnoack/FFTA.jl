@enum Direction FFT_FORWARD=-1 FFT_BACKWARD=1
@enum Pow248 POW2=2 POW4=1 POW8=0
@enum FFTEnum compositeFFT dft pow2FFT pow3FFT pow4FFT pow8FFT

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
Check if `N` is a power of 2, 4, or 8

Returns POW8 (0) if N is a power of 8, POW4 (1) if power of 4, POW2 (2) if power of 2, or nothing otherwise.
"""
function _ispow248(N::Int)
    N < 1 && return nothing
    # Check for powers of 8
    while N & 0b111 == 0
        N >>= 3
    end
    if N == 1
        return POW8
    end
    # Check for powers of 4
    while N & 0b11 == 0
        N >>= 2
    end
    if N == 1
        return POW4
    end
    # Check for powers of 2
    while N & 0b1 == 0
        N >>= 1
    end
    return N == 1 ? POW2 : nothing
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

"""
function CallGraphNode!(nodes::Vector{CallGraphNode{T}}, N::Int, workspace::Vector{Vector{T}}, s_in::Int, s_out::Int)::Int where {T}
    if N == 0
        throw(DimensionMismatch("array has to be non-empty"))
    end
    w = cispi(T(2)/N)
    if iseven(N)
        pow = _ispow248(N)
        if !isnothing(pow)
            if pow == POW8
                # Power of 8: use pow8FFT directly
                push!(workspace, T[])
                push!(nodes, CallGraphNode(0, 0, pow8FFT, N, s_in, s_out, w))
                return 1
            elseif pow == POW4
                # Power of 4: use pow4FFT directly
                push!(workspace, T[])
                push!(nodes, CallGraphNode(0, 0, pow4FFT, N, s_in, s_out, w))
                return 1
            else
                # POW2 - odd power of 2
                if N == 2
                    # Base case: use pow2FFT directly
                    push!(workspace, T[])
                    push!(nodes, CallGraphNode(0, 0, pow2FFT, N, s_in, s_out, w))
                    return 1
                else
                    # Odd power of 2 with N > 2: create composite with pow4FFT + pow2FFT
                    # N = (N/2) * 2, where N/2 is a power of 4
                    N1 = N ÷ 2  # This will be a power of 4
                    N2 = 2       # Remaining power
                    push!(nodes, CallGraphNode(0, 0, dft, N, s_in, s_out, w))
                    sz = length(nodes)
                    push!(workspace, Vector{T}(undef, N))
                    left_len = CallGraphNode!(nodes, N1, workspace, N2, N2*s_out)
                    right_len = CallGraphNode!(nodes, N2, workspace, N1*s_in, 1)
                    nodes[sz] = CallGraphNode(1, 1 + left_len, compositeFFT, N, s_in, s_out, w)
                    return 1 + left_len + right_len
                end
            end
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
        push!(workspace, T[])
        push!(nodes, CallGraphNode(0, 0, dft, N, s_in, s_out, w))
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
    left_len = CallGraphNode!(nodes, N1, workspace, N2, N2*s_out)
    right_len = CallGraphNode!(nodes, N2, workspace, N1*s_in, 1)
    nodes[sz] = CallGraphNode(1, 1 + left_len, compositeFFT, N, s_in, s_out, w)
    return 1 + left_len + right_len
end

"""
$(TYPEDSIGNATURES)
Instantiate a CallGraph from a number `N`

"""
function CallGraph{T}(N::Int) where {T}
    nodes = CallGraphNode{T}[]
    workspace = Vector{Vector{T}}()
    CallGraphNode!(nodes, N, workspace, 1, 1)
    CallGraph(nodes, workspace)
end
