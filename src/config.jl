### possible interpolation modes for EMD.
@enum InterpMode Linear Pchip Cubic Spline

"""
    EMDConfig(Opts...)

EMDConfig sets up the configuration for a decomposition. This
does not include the input signal, but purely the options surrounding
sifting, maximum iterations, the type of interpolation between extrema, etc.
"""
struct EMDConfig{T<:AbstractFloat}
    maxiters::Int
    maxmodes::Int
    ndirs::Int
    sift::Int
    sift_h::Int
    taxis::StepRange
    stop::AbstractVector{T}
    x::AbstractVector{T}
    N::Int
    interp::InterpMode
    function EMDConfig(x::AbstractVector{T}) where T
        N = length(x)
        t = StepRange(1, 1, N)
        new{T}(Int(2_000), typemax(Int), Int(4), 0, 0, t, Array{T}([0.05, 0.5, 0.05]), x, N, Cubic)
    end
end

"""
    (tmin, tmax, zmin, zmax) = boundscheck(opts::EMDConfig, indmin::Array, indmax::Array, nbsym::Int)

Boundary check to define extrema beyond the input signal limits to prevent boundary issues.
Without mirror symmetry, we get ramping at either ends of the IMF signal.
"""
function boundscheck(indmin::AbstractVector, indmax::AbstractVector, opts::EMDConfig, z, nbsym::Int)
    length(indmin) + length(indmax) ≥ 3 || error("EMD: boundscheck: not enough extrema.")
    if (indmax[1] < indmin[1])
        if (opts.x[1] > opts.x[indmin[1]])
            # original EMD-m file used "end" keyword with minimum function??
        end
    end
end

"""
   stop = stopemd(opts::EMDConfig, imf::AbstractVector)

Returns a flag indicating if at least 3 extrema are present to continue 
the decomposition.
"""
function stopemd(otps::EMDConfig, imf::AbstractVector)
    (indmin, indmax) = extr(imf)
    Bool(length(indmin) + length(indmax) < 3)
end


function stopsift(opts::EMDConfig, imf::AbstractVector)
    (mean, npeaks, ndips, amp) = meanamplitude(opts, imf)
    Sx = abs.(mean) ./ amp
    S = mean(Sx)
end

"""
   oind = orthoindex(opts::EMDConfig, imf::AbstractVector)

`orthoindex` computes the index of orthogonality based on the input
signal `x` and the prospective mode function `imf`.
"""
function orthoindex(opts::EMDConfig, imf::AbstractVector)
    n = size(imf, 1)
    s = 0
    for i = 1:n
        for j = 1:n
            if i != j
                # we can improve this with `muladd`...
                s += abs(sum(imf[i,:]) .* conj.(imf[j,:]))/sum(x.^2)
            end
        end
    end
    0.5*s
end
    
    
