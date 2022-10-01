"""
    y = fliplr(x)

alias to flip the input vector `x`.
"""
fliplr(x::AbstractVector) = x[end:-1:1]

"""
    (stop, tmin, tmax, zmin, zmax) = boundarycheck(indmin::Array, indmax::Array, t::Array, z::Array, nbsym::Int)

Boundary check to define extrema beyond the input signal limits to prevent boundary issues.
Without mirror symmetry, we get ramping at either ends of the IMF signal.
"""
function boundarycheck(
    indmin::AbstractVector{T},
    indmax::AbstractVector{T},
    t::AbstractVector{T},
    x::AbstractVector{F},
    z::AbstractVector{F},
    nbsym::Int
) where {T<:Integer,F<:AbstractFloat}
    stop = false
    lx = length(x)
    lsym, rsym = zeros(F, 2)
    if (length(indmin) + length(indmax) < 3)
        tmin, tmax, zmin, zmax = zeros(F, 4)
        stop = true
        return stop, tmin, tmax, zmin, zmax
    end
    if (indmax[1] < indmin[1])
        if (x[1] > x[indmin[1]])
            lmax = fliplr(indmax[2:min(end, nbsym + 1)])
            lmin = fliplr(indmin[1:min(end, nbsym)])
            lsym = indmax[1]
        else
            lmax = fliplr(indmax[1:min(end, nbsym)])
            lmin = vcat(fliplr(indmin[1:min(end, nbsym - 1)]), 1)
            lsym = 1
        end
    else
        if (x[1] < x[indmax[1]])
            lmax = fliplr(indmax[1:min(end, nbsym)])
            lmin = fliplr(indmin[2:min(end, nbsym + 1)])
            lsym = indmin[1]
        else
            lmax = vcat(fliplr(indmax[1:min(end, nbsym - 1)]), 1)
            lmin = fliplr(indmin[1:min(end, nbsym)])
            lsym = 1
        end
    end
    if (indmax[end] < indmin[end])
        if (x[end] < x[indmax[end]])
            rmax = fliplr(indmax[max(end - nbsym + 1, 1):end])
            rmin = fliplr(indmin[max(end - nbsym, 1):end-1])
            rsym = indmin[end]
        else
            rmax = vcat(lx, fliplr(indmax[max(end - nbsym + 2, 1):end]))
            rmin = fliplr(indmin[max(end - nbsym + 1, 1):end])
            rsym = lx
        end
    else
        if (x[end] > x[indmin[end]])
            rmax = fliplr(indmax[max(end - nbsym, 1):end-1])
            rmin = fliplr(indmin[max(end - nbsym + 1, 1):end])
            rsym = indmax[end]
        else
            rmax = fliplr(indmax[max(end - nbsym + 1, 1):end])
            rmin = vcat(lx, fliplr(indmin[max(end - nbsym + 2, 1):end]))
            rsym = lx
        end
    end
    tlmin = @. 2 * t[lsym] - t[lmin]
    tlmax = @. 2 * t[lsym] - t[lmax]
    trmin = @. 2 * t[rsym] - t[rmin]
    trmax = @. 2 * t[rsym] - t[rmax]

    if (tlmin[1] > t[1]) || (tlmax[1] > t[1])
        if (lsym == indmax[1])
            lmax = fliplr(indmax[1:min(end, nbsym)])
        else
            lmin = fliplr(indmin[1:min(end, nbsym)])
        end
        lsym = 1
        tlmin = @. 2 * t[lsym] - t[lmin]
        tlmax = @. 2 * t[lsym] - t[lmax]
    end

    if (trmin[end] < t[lx]) || (trmax[end] < t[lx])
        if (rsym == indmax[end])
            rmax = fliplr(indmax[max(end - nbsym + 1, 1):end])
        else
            rmin = fliplr(indmin[max(end - nbsym + 1, 1):end])
        end
        rsym = copy(lx)
        trmin = @. 2 * t[rsym] - t[rmin]
        trmax = @. 2 * t[rsym] - t[rmax]
    end

    tmin = vcat(tlmin, t[indmin], trmin)
    tmax = vcat(tlmax, t[indmax], trmax)
    zmin = vcat(z[lmin], z[indmin], z[rmin])
    zmax = vcat(z[lmax], z[indmax], z[rmax])

    stop, tmin, tmax, zmin, zmax
end

"""
   stop = stopemd(imf::AbstractVector)

Returns a flag indicating if at least 3 extrema are present to continue 
the decomposition.
"""
function stopemd(imf::AbstractVector{T}) where {T}
    (indmin, indmax) = extrminmax(imf)
    Bool(length(indmin) + length(indmax) < 3)
end

"""
    (stop, envmean, muval) = stopsifting(imf, σ, σ₂, tol)

Default stopping criterion for sifting.
"""
function stopsifting(imf::AbstractVector{T}, σ::T, σ₂::T, tol::T; order::Int=3) where {T<:AbstractFloat}
    (stop, envmean, numextr, numzer, amp) = meanamplitude(imf, order=order)
    if (stop)
        muval = 0
        return stop, envmean, muval
    end
    Sx = @. abs(envmean) / amp
    muval = mean(Sx)
    flag1 = mean(Sx .> σ) > tol
    flag2 = any(Sx .> σ₂)
    flag3 = (numextr > 2)
    flag4 = (abs(numzer - numextr) > 1)
    stop = !((flag1 | flag2) & flag3) && !(flag4)
    stop, envmean, muval
end



"""
   oind = orthoindex(imf::AbstractVector)

`orthoindex` computes the index of orthogonality based on the input
signal `x` and the prospective mode function `imf`.
"""
function orthoindex(x::AbstractVector)
    n = length(x)
    s = 0.0
    for i in eachindex(n)
        for j in eachindex(n)
            if i != j
                # we can improve this with `muladd`...
                s += abs(sum(x[i, :]) .* conj.(x[j, :])) / sum(x .^ 2)
            end
        end
    end
    0.5 * s
end

"""
    (indmin, indmax) = extrminmax(x)
    
`extr` extracts the indices of extrema in value vector `x` over
domain `t`. min/max comparison attempts to match MATLAB behavior.
"""
function extrminmax(x::AbstractVector{T}) where {T<:AbstractFloat}
    dx = diff(x)
    m, n = length(x), length(dx)
    d1 = dx[1:n-1]
    d2 = dx[2:n]
    bad = (dx .== 0)
    indmin = findall((d1 .* d2 .< 0) .& (d1 .< 0)) .+ 1
    indmax = findall((d1 .* d2 .< 0) .& (d1 .> 0)) .+ 1
    if any(bad)
        imax, imin = Int[], Int[]
        dd = diff(vcat(false, bad, false))
        head = findall(dd .== 1)
        tail = findall(dd .== -1)
        if (head[1] == 1)
            if length(head) > 1
                head = head[2:end]
                tail = tail[2:end]
            else
                head, tail = empty(head), empty(tail)
            end
        end
        if length(head) > 0
            if tail[end] == m
                if length(head) > 1
                    head = head[1:end-1]
                    tail = tail[1:end-1]
                else
                    head, tail = empty(head), empty(tail)
                end
            end
        end
        lh = length(head)
        if lh > 0
            for k = 1:lh
                # attempting to match rounding in MATLAB. 0.5 corner-case.
                if dx[head[k]-1] > 0
                    if dx[tail[k]] < 0
                        half = round((head[k] + tail[k]) / 2, RoundNearestTiesUp)
                        append!(imax, half)
                    end
                else
                    if dx[tail[k]] > 0
                        half = round((head[k] + tail[k]) / 2, RoundNearestTiesUp)
                        append!(imin, half)
                    end
                end
            end
        end
        if length(imax) > 0
            sort!(append!(indmax, imax))
        end
        if length(imin) > 0
            sort!(append!(indmin, imin))
        end
    end
    indmin, indmax
end

"""
    (stop, mean, numextr, numzer, amp) = meanamplitude(x::AbstractVector, order::Int=3)

Computes the mean of the envelopes, the mode amplitude estimate, and the number
of extrema, including zeros. String should indicate "cubic", "linear", etc.
"""

function meanamplitude(x::AbstractVector{T}; order::Int=3) where {T<:AbstractFloat}
    N = length(x)
    envmin, envmax, envmean, amp = zeros(T, N), zeros(T, N), zeros(T, N), zeros(T, N)
    amp = zero(T)
    stop = false
    t = collect(1:N)
    (indmin, indmax) = extrminmax(x)
    indzer = extrzeros(x)
    numextr = length(indmin) + length(indmax)
    numzer = length(indzer)
    (stop, tmin, tmax, mmin, mmax) = boundarycheck(
        indmin, indmax, t, x, x, 2)
    if (stop)
        return stop, mean, numextr, numzer, amp
    end

    # construct the interpolant and then pass the x-axis. Corner check for end-of-knot condition.
    # we can use parallel programming for interpolating the min/max envelopes.
    @sync begin
        @async envmin = ifelse((length(mmin) == 3), parabolaspline(tmin, mmin, t), Spline1D(tmin, mmin; k=order)(t))
        @async envmax = ifelse((length(mmax) == 3), parabolaspline(tmax, mmax, t), Spline1D(tmax, mmax; k=order)(t))
    end

    # get the μ-envelope from the extrema splines.
    envmean = @. (envmin + envmax) / 2

    # in MATLAB, this was mean(abs(envmax-envmin),1)/2, but since dim=1,
    # this is effectively a meaningless call since they're scalars.
    amp = @. abs(envmax - envmin) / 2
    stop, envmean, numextr, numzer, amp
end


"""
    Vq = parabolaspline(x::AbstractArray, v::AbstractArray, Xq::AbstractArray)

Computes cubic spline with parabolic data and end-of-knot condition.
See De Boor, C. (1978). A practical guide to splines (Vol. 27, p. 325). 
New York: springer-verlag.
"""
function parabolaspline(X::AbstractVector, V::AbstractVector, Xq::AbstractVector{T}) where {T}
    M = length(Xq)
    Vq = zeros(T, M)
    dx, dv = diff(X), diff(V)
    dx2 = X[3] - X[1] # skipping second "tau" point
    dvdx = dv ./ dx
    coeffs = zeros(3) # ax² + bx + c
    coeffs[1] = (diff(dvdx)/dx2)[1]
    coeffs[2] = dvdx[1] - coeffs[1] * dx[1]
    coeffs[3] = V[1]
    Vq = fill(coeffs[1], M)
    xs = collect(Xq) .- fill(X[1], M)
    for i = 2:3
        Vq = Vq .* xs + fill(coeffs[i], M)
    end
    Vq
end

"""
    indzers = extrzeros(x)

Extract indices of zero-crossing extrema.
"""
function extrzeros(x::AbstractVector)
    x1 = x[1:end-1]
    x2 = x[2:end]
    indzer = findall(x1 .* x2 .< 0)
    iz = (x .== 0)
    if any(iz)
        zer = findall(x .== 0)
        if any(diff(iz) .== 1)
            dz = diff(vcat(false, iz, false))
            headz = findall(dz .== 1)
            tailz = findall(dz .== -1) .- 1
            indz = @. (headz + tailz) / 2
            indz = round.(Int, indz, RoundNearestTiesUp)
        else
            indz = iz
        end
        sort!(append!(indzer, indz))
    end
    indzer
end
