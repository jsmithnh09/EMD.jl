using Pkg
Pkg.activate(".")
using EMD


# first sifting iteration variables stored in MATLAB, comparing Spline1D operation.
fpath = joinpath(@__DIR__, "signals")

x = Array{Float64, 1}(undef, 4096)
read!(joinpath(fpath, "doppler.bin"), x)

mat_zmin = Vector{Float64}(undef, 24)
mat_zmax = Vector{Float64}(undef, 24)
read!(joinpath(fpath, "zmin.bin"), mat_zmin)
read!(joinpath(fpath, "zmax.bin"), mat_zmax)

mat_imf = Vector{Float64}(undef, 4096*4)
read!(joinpath(fpath, "imf4x4096.bin"), mat_imf)
mat_imf = reshape(mat_imf, (4096, 4))

lx = length(x)
t = collect(1:lx)
fliplr(x) = x[end:-1:1]
nbsym = 2
z = copy(x)

# manually confirmed this is correct.
(indmin, indmax) = EMD.extrminmax(x)


# first, check we're extracting the correct lmin/lmax/rmin/rmax indices.
if (indmax[1] < indmin[1])
    if (x[1] > x[indmin[1]])
        lmax = fliplr(indmax[2:min(end, nbsym+1)])
        lmin = fliplr(indmin[1:min(end, nbsym)])
        lsym = copy(indmax[1])
    else
        lmax = fliplr(indmax[1:min(end, nbsym)])
        lmin = vcat(fliplr(indmin[1:min(end, nbsym-1)]), 1)
        lsym = 1
    end
else
    if (x[1] < x[indmax[1]])
        lmax = fliplr(indmax[1:min(end, nbsym)])
        lmin = fliplr(indmin[2:min(end, nbsym+1)])
        lsym = copy(indmin[1])
    else
        lmax = vcat(fliplr(indmax[1:min(end, nbsym-1)]), 1)
        lmin = fliplr(indmin[1:min(end, nbsym)])
        lsym = 1
    end
end
if (indmax[end] < indmin[end])
    if (x[end] < x[indmax[end]])
        rmax = fliplr(indmax[max(end-nbsym+1, 1):end])
        rmin = fliplr(indmin[max(end-nbsym, 1):end-1])
        rsym = copy(indmin[end])
    else
        rmax = vcat(lx, fliplr(indmax[max(end-nbsym+2, 1):end]))
        rmin = fliplr(indmin[max(end-nbsym+1, 1):end])
        rsym = copy(lx)
    end
else
    if (x[end] > x[indmin[end]])
        rmax = fliplr(indmax[max(end-nbsym, 1):end-1])
        rmin = fliplr(indmin[max(end-nbsym+1, 1):end])
        rsym = copy(indmax[end])
    else
        rmax = fliplr(indmax[max(end-nbsym+1, 1):end])
        rmin = vcat(lx, fliplr(indmin[max(end-nbsym+2, 1):end]))
        rsym = copy(lx)
    end
end

zlmax = z[lmax]
zlmin = z[lmin]
zrmax = z[rmax]
zrmin = z[rmin]

zmin = vcat(zlmin, z[indmin], zrmin)
zmax = vcat(zlmax, z[indmax], zrmax)