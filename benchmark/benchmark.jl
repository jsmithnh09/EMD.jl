using BenchmarkTools, Wavelets, EMD
using HDF5

include("noise.jl")

mutable struct EMDResult
    snr::Int
    size::Int
    time::Float64
    signal::String
end

## Setting up test conditions.
# snr = -10:2:10 # ignore the SNR dimension of testing for now.
N = 8:2:16
precompile(emd, (Vector{Float64},))
snames = (:Doppler, :HeaviSine, :Blocks, :Bumps)

#=
# run first-time JIT (this doesn't need to be done if using BenchmarkTools.jl since the function is called several times.)
# _ = emd(zeros(4097))
=#

## pre-alloc the vector Array based on the  parameters above
# results = Vector{EMDResult}(undef, length(N) * length(snr) * length(snames))
h5open(joinpath(pwd(), "emd_benchmark.hdf5"), "w") do fid
    g = create_group(fid, "env")
    write(g, "noise_axis", collect(Int64, N))
    write(g, "signals", ["$name" for name in snames])

    for signal in snames
        sg = create_group(fid, "$(signal)_results")
        t_res = zeros((length(N),))
        ind = 1
        for n in N
            X = testfunction(2^n, "$signal")
            # acn!(X, s, "pink") # ignoring SNR - only looking at timing with increasing N.
            b = @benchmark emd(X)
            t_res[ind] = minimum(b).time / 1e9 # BenchmarkTools stores in nanoseconds. 
            ind += 1
        end
        write(sg, "time", t_res)
    end
end

