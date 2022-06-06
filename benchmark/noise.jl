using FFTW
using Statistics: mean, std

export acn, acn!

const ALPHA_FACTORS = Dict{String, Int8}(
    "white" => Int8(0),
    "pink" => Int8(-1),
    "violet" => Int8(2),
    "red" => Int8(-2),
    "blue" => Int8(1)
)

acn!(x::AbstractVector, dBsnr::Int, color::String="white") = acn!(x, Float64(dBsnr), color)
function acn!(x::AbstractVector, dBsnr::AbstractFloat, color::String="white")
    linsnr = 10^(dBsnr / 10)
    α = ALPHA_FACTORS[color]
    # convert alpha-factor from PSD to Amplitude Density.
    α /= 2
    N = length(x)
    # generating white noise; don't manipulate spectra if white noise.
    w = randn(N)
    # manipulate spectral roll-off for left-side of spectrum.
    if (color != "white")
        M = ceil(Int, (N/2)+1)
        W = fft(w, (1,))
        fv = collect(Int, 1:M)
        W = W[fv]
        W = W .* (fv.^α)
        # even or odd, in/exclude Nyquist when mirroring.
        if (rem(N,2) == 1)
            W = vcat(W, conj.(W[end:-1:2]))
        else
            W = vcat(W, conj.(W[end-1:-1:2]))
        end
        n = real.(ifft(W, (1,)))
    else
        n = w
    end
    # zero mean, unity variance noise
    n .-= mean(n)
    n ./= std(n)
    # scale the noise based on SNR
    Ps, Pn = sum(x.^2), sum(n.^2)
    β = √(Ps / (linsnr * Pn))
    x .+= (β .* n)
    x
end

function acn(x::AbstractVector, dBsnr::AbstractFloat, color::String="white")
    y = copy(x)
    acn!(y, dBsnr, color)
    y
end

