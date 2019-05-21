using Interpolations
using Dierckx
using Statistics: mean

"""
    fliplr(x)

array returned with order of elements flipped from left to right.
"""
function fliplr(x)
    input = copy(x)
    output = input[end:-1:1]
    return output
end # fliplr

#-------------------------------------------------------------------------------

"""
    interp1(xpt, ypt, x; method, extrapValue)

MATLAB-similar syntax for interpolation. This matches an example provided
[here](https://discourse.julialang.org/t/julia-version-of-matlab-function-interp1/8540).
Alternatively, `Spline1D(x, v; k=3)` can be used for cubic splining. As a single line,
this would be `Spline1D(x, v; k=3)(xq)`.
"""
function interp1(xpt, ypt, x; method="cubic", extrapValue=nothing)

    if extrapValue == nothing
        y = zeros(Float64, length(x))            # pre-defining type as Float.
        idx = trues(length(x))                   # BitArray, all set to true
    else
        y = extrapValue*ones(length(x))
        idx = (x .>= xpt[1]) .& (x .<= xpt[end]) # index within allowable range
    end

    if method == "linear"
        intf = interpolate((xpt,), ypt, Gridded(Linear()))
        y[idx] = intf[x[idx]]

    elseif method == "cubic"
        # interpolations object.
        itp = interpolate(ypt, BSpline(Cubic(Natural())), OnGrid())
        intf = scale(itp, xpt)

        # populate with interpolated objects
        y[idx] = [intf[xi] for xi in x[idx]]
    end

    return y
end
