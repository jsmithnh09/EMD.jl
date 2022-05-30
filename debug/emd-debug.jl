using Pkg
using Plots
Pkg.activate(".")
using EMD

fpath = joinpath(@__DIR__, "signals")
X = Vector{Float64}(undef, 4096)
read!(joinpath(fpath, "doppler.bin"), X)

mat_imf = Vector{Float64}(undef, 4096*4)
read!(joinpath(fpath, "imf4096x4.bin"), mat_imf)
mat_imf = reshape(mat_imf, (4096, 4))

# perform the default EMD algorithm, no modifiers.
imf = EMD.emd(X)
