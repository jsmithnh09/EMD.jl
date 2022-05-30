## REPL-flavored EMD implementation.

using Pkg
Pkg.activate(".")
cd("src")
using EMD
using DelimitedFiles
using Wavelets: testfunction

T = Float64

# pull in the components. 2^12 should yield 3 IMFs.
x = vec(readdlm("heavysine12.txt", T))
x = testfunction(2^14, "Doppler")

maxmodes = 0
maxiters = Int(2_000)


### Main EMD loop
k, curiter = 1, 0
r = copy(x)

imf = Vector{Vector{T}}()

while(!EMD.stopemd(r) && (maxmodes == 0 || k < maxmodes + 1))
    m = copy(r)
    (stopsift, μenv, _) = EMD.stopsifting(m, 0.05, 0.5, 0.05)

    # sift loop
    while ((!stopsift) && (curiter < maxiters))
        m -= μenv
        (stopsift, μenv, _) = EMD.stopsifting(m, 0.05, 0.5, 0.05)
        curiter += 1
    end
    push!(imf, m)
    k += 1
    r -= m
    curiter = 0
end
if (any(r != 0))
    push!(imf, r)
end
