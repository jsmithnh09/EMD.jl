using EMD, Test

@testset "EMD utilities" begin
    include("util.jl")
end

# deconstruct the signal into IMFs, then sum to prove reconstruction.
@testset "Partial Reconstruction" begin
    npts = 2^14
    
    # just using HeaviSine as a test input signal from
    #   DL Donoho and IM Johnstone. Ideal spatial separation by wavelet shrinkage.
    #   Biometrika, vol. 81, pp.425-455, 1994.
    
    x = Vector{Float64}(undef, npts)
    taxis = 0:1/npts:1-eps()
    i = 1
    for t in eachindex(taxis)
        x[i] = 4*sin(4π*t) - sign(t-0.3) - sign(0.72-t)
        i += 1
    end
    imf = emd(x)
    @test isapprox(sum(imf), x, rtol=1e-8)
end

# setting up conditionals prior to testing full algorithm.
fpath = joinpath(@__DIR__, "m-data")
X = Vector{Float64}(undef, 4096)
read!(joinpath(fpath, "doppler.bin"), X)
mat_imf = Vector{Float64}(undef, 4 * 4096) # encoded 4-IMF into column major, emd_rilling gives back row-major.
read!(joinpath(fpath, "imf4096x4.bin"), mat_imf)
mat_imf = reshape(mat_imf, (4096, 4))

# determine if the sifted IMF's match the MATLAB result.
@testset "Primary EMD-algorithm" begin
    imf = EMD.emd(X)
    for ind in eachindex(imf)
        @test isapprox(imf[ind], mat_imf[:, ind], rtol=1e-8) # √ϵ
    end
end