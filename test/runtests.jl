using DelimitedFiles

# % [t, X] with X = wnoise(4,12)
# fwrite(fid, '%2.1f %12.8f \n', A);

# return carriage or EOL creates the row entry.
X = readdlm("doppler.txt")

# perform EMD and reconstruct the signal to variable Y.
for i = 1:size(X, 1)
    X[i, 1] ≈ Y[i, 2] atol = 1e-5
end
