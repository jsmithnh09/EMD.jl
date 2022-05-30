using EMD, Test
using Wavelets

@testset "EMD utilities" begin
    include("util.jl")
end

# deconstruct the signal into IMFs, then sum to prove reconstruction.
@testset "Partial Reconstruction" begin
    for sig in (:Doppler, :HeaviSine, :Blocks, :Bumps)
        @testset "$(sig) Signal" begin
            x = testfunction(2^14, "$(sig)")
            imf = emd(x)
            @test isapprox(sum(imf), x, rtol=1e-8)
        end
    end
end

# setting up conditionals prior to testing full algorithm.
fpath = joinpath(@__DIR__, "m-data")
X = Vector{Float64}(undef, 4096)
read!(joinpath(fpath, "doppler.bin"), X)
mat_imf = Vector{Float64}(undef, 4*4096) # encoded 4-IMF into column major, emd_rilling gives back row-major.
read!(joinpath(fpath, "imf4096x4.bin"), mat_imf)
mat_imf = reshape(mat_imf, (4096, 4))

# determine if the sifted IMF's match the MATLAB result.
@testset "Primary EMD-algorithm" begin
    imf = EMD.emd(X)
    for ind = 1:length(imf)
        @test isapprox(imf[ind], mat_imf[:,ind], rtol=1e-8) # √ϵ
    end    
end