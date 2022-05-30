# EMD.jl
This project aims to implement [Empirical Mode Decomposition](https://en.wikipedia.org/wiki/Hilbert-Huang_transform), generating Intrinsic Mode Functions (IMFs) [1]. This project adapts the original C and MATLAB script provided by Rilling _et al_ [2] into the Julia Language, where sifting stop criterion are defined and window boundary effects are properly addressed.

_Contributors_: Jordan R. Smith

[1] N.E. Huang et. al., "The empirical mode decomposition and the Hilbert spectrum for non-linear stationary time series analysis", Proc. Royal Soc. London A, Vol. 454, pp. 903-995, 1998

[2] G. Rilling, P. Flandrin, and P. Goncalves "On Empirical Mode Decomposition and its algorithms", IEEE-EURASIP Workshop on Nonlinear Signal and Image Processing NSIP-03, Grado (I), June 2003

## TODO
- Re-use the time-axis since its always the same `UnitRange` type
- Pre-allocate and re-use the envelope min/max buffers since they'll be the same length
- Pass the Config type to the utility functions to save on any memcopies
