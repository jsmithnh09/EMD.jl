using BenchmarkTools, Wavelets, EMD
using SpecialFunctions: erfc

include("noise.jl")

mutable struct EMDResult 
    snr::Int
    size::Int
    time::Float64
    signal::String
end

## Setting up test conditions.
N = 8:2:16
snr = -10:2:10
precompile(emd, (Vector{Float64},))
snames = (:Doppler, :HeaviSine, :Blocks, :Bumps)

#=
# run first-time JIT (this doesn't need to be done if using BenchmarkTools.jl since the function is called several times.)
# _ = emd(zeros(4097))
=#

## pre-alloc the vector Array based on the  parameters above
results = Vector{EMDResult}(undef, length(N)*length(snr)*length(snames))
curind = 1
for signal in snames
    for n in N
        X = testfunction(2^n, "$signal")
        for s in snr
            acn!(X, s, "pink")
            b = @benchmark emd(X)
            t = minimum(b).time / 1e9 # BenchmarkTools stores in nanoseconds.
            results[curind] = EMDResult(s, n, t, "$signal")
            curind += 1