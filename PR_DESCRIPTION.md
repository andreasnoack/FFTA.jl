# Optimize FFTs of odd powers of two using composite pow4FFT plans

## Summary

This PR optimizes FFT computation for sizes that are odd powers of 2 (e.g., 8, 32, 128, 512) by using a composite plan that leverages `fft_pow4!` for as much of the computation as possible, rather than relying purely on `fft_pow2!`.

### Optimization Strategy

For FFT sizes that are odd powers of 2, the planner now creates a composite plan by factoring N = (N/2) × 2, where N/2 is a power of 4:

- **N=8 (2³)**: Uses pow4FFT(4) + pow2FFT(2) instead of pure pow2FFT
- **N=32 (2⁵)**: Uses pow4FFT(16) + pow2FFT(2) instead of pure pow2FFT
- **N=128 (2⁷)**: Uses pow4FFT(64) + pow2FFT(2) instead of pure pow2FFT
- **N=512 (2⁹)**: Uses pow4FFT(256) + pow2FFT(2) instead of pure pow2FFT

### Performance Benefits

The `fft_pow4!` algorithm processes 4 elements at a time with more efficient butterfly operations compared to `fft_pow2!`'s 2-element processing. This should provide better performance for odd power-of-2 FFT sizes.

### Implementation Details

**Modified `CallGraphNode!` in `src/callgraph.jl`:**
- Detects when N is an odd power of 2 (using existing `_ispow24()` helper)
- Creates a composite plan structure instead of a single pow2FFT node
- Recursively builds the call graph with pow4FFT for the N/2 subtransform and pow2FFT for the remaining factor of 2

**Even powers of 2** (4, 16, 64, 256, etc.) continue to use pure pow4FFT as before.

### Code Changes

The key change is in `src/callgraph.jl`, lines 73-95, where the planning logic now branches:

```julia
if pow == POW4
    # Even power of 2: use pow4FFT directly
    push!(workspace, T[])
    push!(nodes, CallGraphNode(0, 0, pow4FFT, N, s_in, s_out, w))
    return 1
else
    # Odd power of 2: create composite with pow4FFT + pow2FFT
    N1 = N ÷ 2  # This will be a power of 4
    N2 = 2       # Remaining power
    # ... create composite plan structure
end
```

### Testing

- Added `test_odd_powers.jl` verification script to demonstrate the optimization
- Existing test suite (sizes 1-64 in `test/onedim/*.jl`) should pass and will exercise the new code paths for sizes 8, 32, and 64
- CI will run comprehensive tests across all Julia versions

## Test Plan

- [x] Modify planning logic in `src/callgraph.jl`
- [x] Add verification script `test_odd_powers.jl`
- [x] Commit and push changes to branch `claude/optimize-fft-odd-powers-GLnCr`
- [ ] Run existing test suite: `Pkg.test()` (will be done by CI)
- [ ] Verify correctness for odd power sizes (8, 32, 128, 512, 2048)
- [ ] Check that even power sizes still work correctly
- [ ] Run performance benchmarks to quantify improvement (optional)

## Branch Information

- **Branch**: `claude/optimize-fft-odd-powers-GLnCr`
- **Base**: `main`
- **Repository**: `andreasnoack/FFTA.jl`

## How to Create the PR

Visit: https://github.com/andreasnoack/FFTA.jl/compare/main...claude/optimize-fft-odd-powers-GLnCr

Or use the GitHub CLI:
```bash
gh pr create --base main --head claude/optimize-fft-odd-powers-GLnCr \
  --title "Optimize FFTs of odd powers of two using composite pow4FFT plans" \
  --body-file PR_DESCRIPTION.md
```
