using BenchmarkTools
using EMD
using Wavelets

include("noise.jl")

mutable struct EMDResult 
    snr::Int
    size::Int
    result::BenchmarkTools.Trial
end


N = 2:2:18
snr = -14:2:12
precompile(emd, (Vector{Float64},))

for signal in (:Doppler, :HeaviSine, :Blocks, :Bumps)
    for n in N
        X = testfunction(2^n, "$(signal)")
        for s in snr
            
        