"""
Common benchmark definitions shared between FFTA and FFTW benchmarks
"""

using Primes

# Cumulative products of 3,4,5,5,7,11
const CUMULATIVE_PRODUCTS = cumprod([3, 4, 5, 5, 7, 11])

# Select 20 primes with logarithmic spacing
const SELECTED_PRIMES = begin
    # Generate logarithmically spaced target points
    all_primes = primes(20000)
    log_points = exp.(range(log(2), log(20000), length=20))

    # Find nearest prime to each logarithmic point
    selected = Int[]
    for target in log_points
        # Find the prime closest to this target
        idx = argmin(abs.(all_primes .- target))
        candidate = all_primes[idx]
        # Avoid duplicates
        if !(candidate in selected)
            push!(selected, candidate)
        end
    end

    sort!(selected)
end

# Odd powers of 2 (up to 2^19 = 524288)
const ODD_POWERS_OF_2 = [2^i for i in 1:2:19]

# Even powers of 2 (up to 2^20 = 1048576)
const EVEN_POWERS_OF_2 = [2^i for i in 2:2:20]

# Powers of 3
const POWERS_OF_3 = [3^i for i in 1:9]

"""
    create_size_categories()

Create a dictionary mapping array sizes to their categories.
"""
function create_size_categories()
    categories = Dict{Int, String}()

    # Categorize odd powers of 2
    for n in ODD_POWERS_OF_2
        categories[n] = "odd_power_of_2"
    end

    # Categorize even powers of 2
    for n in EVEN_POWERS_OF_2
        categories[n] = "even_power_of_2"
    end

    # Categorize powers of 3
    for n in POWERS_OF_3
        categories[n] = "power_of_3"
    end

    # Categorize composite numbers
    for n in CUMULATIVE_PRODUCTS
        categories[n] = "composite"
    end

    # Categorize primes
    for n in SELECTED_PRIMES
        categories[n] = "prime"
    end

    return categories
end

"""
    get_all_sizes()

Get all benchmark sizes in a flat array.
"""
function get_all_sizes()
    return vcat(ODD_POWERS_OF_2, EVEN_POWERS_OF_2, POWERS_OF_3,
                CUMULATIVE_PRODUCTS, SELECTED_PRIMES)
end

# ===== 2D Problem Sizes =====

"""
2D sizes for benchmarking: (rows, cols) tuples
Organized by category to test different combinations
"""

# Square arrays from powers of 2 (limited to reasonable sizes for 2D)
const SQUARE_POWERS_OF_2_2D = [(2^i, 2^i) for i in 4:10]  # 16x16 to 1024x1024

# Square arrays from powers of 3
const SQUARE_POWERS_OF_3_2D = [(3^i, 3^i) for i in 2:6]  # 9x9 to 729x729

# Rectangular arrays: power of 2 x power of 2 (different dimensions)
const RECT_POWER2_2D = [
    (64, 128), (128, 64), (256, 512), (512, 256),
    (128, 1024), (1024, 128)
]

# Mixed rectangular: power of 2 x power of 3
const MIXED_POWER2_POWER3_2D = [
    (64, 81), (81, 64), (128, 243), (243, 128),
    (256, 243), (243, 256)
]

# Mixed: power of 2 x prime
const MIXED_POWER2_PRIME_2D = [
    (64, 127), (127, 64), (128, 251), (251, 128),
    (256, 509), (509, 256)
]

# Mixed: prime x prime (smaller sizes for reasonable performance)
const PRIME_PRIME_2D = [
    (31, 37), (53, 61), (97, 103),
    (127, 131), (251, 257)
]

# Composite sizes
const COMPOSITE_2D = [
    (60, 60), (120, 60), (60, 120),  # 3*4*5
    (84, 84), (168, 84), (84, 168),   # 3*4*7
    (100, 100), (200, 100), (100, 200) # Powers of 10
]

"""
    create_2d_size_categories()

Create a dictionary mapping 2D array sizes to their categories.
"""
function create_2d_size_categories()
    categories = Dict{Tuple{Int,Int}, String}()

    # Categorize square powers of 2
    for size in SQUARE_POWERS_OF_2_2D
        categories[size] = "square_power_of_2"
    end

    # Categorize square powers of 3
    for size in SQUARE_POWERS_OF_3_2D
        categories[size] = "square_power_of_3"
    end

    # Categorize rectangular power of 2
    for size in RECT_POWER2_2D
        categories[size] = "rect_power_of_2"
    end

    # Categorize mixed power 2 x power 3
    for size in MIXED_POWER2_POWER3_2D
        categories[size] = "mixed_power2_power3"
    end

    # Categorize mixed power 2 x prime
    for size in MIXED_POWER2_PRIME_2D
        categories[size] = "mixed_power2_prime"
    end

    # Categorize prime x prime
    for size in PRIME_PRIME_2D
        categories[size] = "prime_prime"
    end

    # Categorize composite sizes
    for size in COMPOSITE_2D
        categories[size] = "composite_2d"
    end

    return categories
end

"""
    get_all_2d_sizes()

Get all 2D benchmark sizes in a flat array of tuples.
"""
function get_all_2d_sizes()
    return vcat(SQUARE_POWERS_OF_2_2D, SQUARE_POWERS_OF_3_2D,
                RECT_POWER2_2D, MIXED_POWER2_POWER3_2D,
                MIXED_POWER2_PRIME_2D, PRIME_PRIME_2D,
                COMPOSITE_2D)
end
