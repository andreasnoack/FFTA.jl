#!/usr/bin/env julia

"""
Benchmark script for FFTW.jl
This script runs in an isolated environment with only FFTW.jl loaded.
"""

# Ensure packages are installed
import Pkg
Pkg.instantiate()

using FFTW
using BenchmarkTools
using JSON

# Load shared benchmark definitions from parent project
# (common_defs.jl uses Primes which is in the parent benchmark project)
const SIZE_CATEGORIES, ALL_SIZES, SIZE_CATEGORIES_2D, ALL_SIZES_2D = let
    parent_project = joinpath(@__DIR__, "..")
    Pkg.activate(parent_project)
    Pkg.instantiate()
    include(joinpath(parent_project, "common_defs.jl"))
    # Restore this project
    Pkg.activate(@__DIR__)
    (create_size_categories(), sort(get_all_sizes()),
     create_2d_size_categories(), get_all_2d_sizes())
end

# Number of samples for each benchmark
const SAMPLES = 100
const EVALS = 10

function benchmark_fftw()
    results = Dict{String, Any}()
    results["package"] = "FFTW"
    results["data"] = []
    results["categories"] = SIZE_CATEGORIES

    println("Benchmarking FFTW.jl...")
    println("=" ^ 50)

    for n in ALL_SIZES
        category = SIZE_CATEGORIES[n]
        println("Testing array size: $n (category: $category)")

        # Benchmark complex FFT
        x = randn(ComplexF64, n)
        trial = @benchmark fft($x) samples=SAMPLES evals=EVALS

        median_time = median(trial).time * 1e-9  # Convert to seconds
        runtime_per_element = median_time / n

        push!(results["data"], Dict(
            "size" => n,
            "category" => category,
            "median_time" => median_time,
            "runtime_per_element" => runtime_per_element,
            "mean_time" => mean(trial).time * 1e-9,
            "min_time" => minimum(trial).time * 1e-9,
            "max_time" => maximum(trial).time * 1e-9
        ))

        println("  Median time: $(median_time * 1e6) μs")
        println("  Time per element: $(runtime_per_element * 1e9) ns")
    end

    # Save results to JSON
    output_file = joinpath(@__DIR__, "..", "results_fftw.json")
    open(output_file, "w") do io
        JSON.print(io, results, 2)
    end

    println("\nResults saved to: $output_file")
    return results
end

function benchmark_fftw_2d()
    results = Dict{String, Any}()
    results["package"] = "FFTW"
    results["dimension"] = "2D"
    results["data"] = []

    # Convert tuple keys to strings for JSON serialization
    categories_serializable = Dict{String, String}()
    for (size_tuple, category) in SIZE_CATEGORIES_2D
        categories_serializable["$(size_tuple[1])x$(size_tuple[2])"] = category
    end
    results["categories"] = categories_serializable

    println("\nBenchmarking FFTW.jl (2D FFTs)...")
    println("=" ^ 50)

    for (rows, cols) in ALL_SIZES_2D
        category = SIZE_CATEGORIES_2D[(rows, cols)]
        total_elements = rows * cols
        println("Testing array size: $(rows)x$(cols) (category: $category, total: $total_elements elements)")

        # Benchmark complex 2D FFT
        x = randn(ComplexF64, rows, cols)
        trial = @benchmark fft($x) samples=SAMPLES evals=EVALS

        median_time = median(trial).time * 1e-9  # Convert to seconds
        runtime_per_element = median_time / total_elements

        push!(results["data"], Dict(
            "rows" => rows,
            "cols" => cols,
            "size_string" => "$(rows)x$(cols)",
            "total_elements" => total_elements,
            "category" => category,
            "median_time" => median_time,
            "runtime_per_element" => runtime_per_element,
            "mean_time" => mean(trial).time * 1e-9,
            "min_time" => minimum(trial).time * 1e-9,
            "max_time" => maximum(trial).time * 1e-9
        ))

        println("  Median time: $(median_time * 1e6) μs")
        println("  Time per element: $(runtime_per_element * 1e9) ns")
    end

    # Save results to JSON
    output_file = joinpath(@__DIR__, "..", "results_fftw_2d.json")
    open(output_file, "w") do io
        JSON.print(io, results, 2)
    end

    println("\nResults saved to: $output_file")
    return results
end

# Run benchmarks
benchmark_fftw()
benchmark_fftw_2d()
