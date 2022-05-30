"""
Empirical Mode Decomposition Package.
"""
module EMD


using Statistics: mean
using Dierckx: Spline1D

export emd

## include environment setup.
include("config.jl")
include("util.jl")

"""
    imf = emd(x; kwargs...)

Given an input vector, Empirical Mode Decomposition is performed. Currently,
only cubic spline interpolation with default stoppage criterion is in place. The return
type `imf` will be a `Vector{Vector{Float64}}`.
"""
function emd(x::AbstractVector{T}; kwargs...) where {T <: AbstractFloat}
    
    # create the configuration object.
    cfg = EMDConfig(x, kwargs...)
    k, curiter = 1, 0
    r = copy(x)
    imf = Vector{Vector{T}}()
    Xtol = √eps(eltype(x)) * maximum(abs, x)


    while (!(stopemd(r)) && (k < cfg.maxmodes+1 || cfg.maxmodes == 0))
        m = copy(r)
        (stopsift, μenv, _) = stopsifting(m, cfg.stop[1], cfg.stop[2], cfg.stop[3], order=cfg.interp)

        # check for noise floor extrema prior to sifting loop.
        if (maximum(abs, m) < Xtol)
            break
        end

        # sift loop
        while ((!stopsift) && (curiter < cfg.maxiters))
            m -= μenv
            (stopsift, μenv, _) = stopsifting(m, cfg.stop[1], cfg.stop[2], cfg.stop[3], order=cfg.interp)
            curiter += 1
        end
        push!(imf, m) # append IMF
        k += 1        # increment the mode iteration.
        r -= m        # extract mode from input.
        curiter = 0   # reset the sifting iteration.
    end
    if (any(r != 0))
        push!(imf, r) # append residual if any fluctuations are still present.
    end
    imf
end # emd

end # module
