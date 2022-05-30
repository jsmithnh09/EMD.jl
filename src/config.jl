"""
    EMDConfig(Opts...)
EMDConfig sets up the configuration for a decomposition. This
does not include the input signal, but purely the options surrounding
sifting, maximum iterations, the type of interpolation between extrema, etc.
"""

const InterpMode = Dict{Int64, String}(
    1 => "linear",
    2 => "quadratic",
    3 => "cubic"
    )

struct EMDConfig
    input::Vector
    maxiters::Int
    maxmodes::Int
    ndirs::Int
    stop::Vector
    interp::Int
    function EMDConfig(
        x::Vector;
        maxiters::Int=2_000, 
        maxmodes::Int=0, 
        ndirs::Int=4, 
        stop::Vector=[0.05, 0.5, 0.05], 
        interp::Int=3
        )
        if length(stop) != 3
            throw(DimensionMismatch("stop vector must be 3-elements, instead got $(size(stop))"))
        elseif (interp < 1 || interp > 3)
            throw(ArgumentError("interpolation order must be in 0 < k ≤ 3."))
        end
        return new(x, maxiters, maxmodes, ndirs, stop, interp)
    end
end

### Show details about the config struct.
function Base.show(io::IO, cfg::EMDConfig)
    mmstr = (cfg.maxmodes == typemax(Int)) ? "∞" : "$(cfg.maxmodes)"
    print(io, "EMD Config:\n") 
    print(io, "maxmodes = ", mmstr, "\n")
    print(io, "stop     = ", cfg.stop, "\n") 
    print(io, "interp   = ", InterpMode[cfg.interp], "\n")
end