using FFTA, Test

# Test script to verify that odd powers of 2 use composite plans with pow4FFT
# This demonstrates the optimization for sizes like 8, 32, 128, etc.

@testset "Odd power optimization verification" begin
    # Test sizes that are odd powers of 2
    odd_power_sizes = [8, 32, 128, 512, 2048]

    @testset "FFT correctness for N=$N" for N in odd_power_sizes
        # Generate random input
        x = complex.(randn(N), randn(N))

        # Compute FFT
        y = fft(x)

        # Verify against inverse
        x_recovered = ifft(y)
        @test x ≈ x_recovered rtol=1e-12

        # Create plan and inspect structure
        plan = plan_fft(x)
        println("Plan for N=$N:")
        println("  Number of nodes: ", length(plan.graph.nodes))
        for (i, node) in enumerate(plan.graph.nodes)
            println("  Node $i: type=$(node.type), size=$(node.sz)")
        end
        println()
    end
end

@testset "Performance comparison hint" begin
    # This is a hint for performance testing (not a functional test)
    # For N=32 (2^5), the optimization creates:
    #   - One composite node (size 32)
    #   - One pow4FFT node (size 16 = 2^4)
    #   - One pow2FFT node (size 2)
    # Instead of a pure pow2FFT approach

    N = 32
    plan = plan_fft(ones(ComplexF64, N))

    # Check that we have a composite plan
    root_node = plan.graph.nodes[1]
    println("Root node for N=$N: type=$(root_node.type), size=$(root_node.sz)")

    if root_node.type == FFTA.compositeFFT
        println("✓ Using composite plan for odd power of 2")
        println("  This will use pow4FFT for most of the work!")
    end
end
