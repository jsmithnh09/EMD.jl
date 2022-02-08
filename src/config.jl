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
    interp::InterpMode
    function EMDConfig(x::AbstractVector{T}) where T <: AbstractFloat
        new{T}(Int(2_000), typemax(Int), 4, 0, 0, StepRange(1, 1, length(x)), Vector{T}([0.05, 0.5, 0.5]), Cubic)
    end
end