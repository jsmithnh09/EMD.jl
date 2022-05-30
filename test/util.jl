using Test, EMD
in = [-4., 3., -2., 0., -4., 1., -3., 2., 2., 3.]
t = collect(1:length(in))

# test that the extrema extraction. using randi([-5, 5], 1, 10)
(minima, maxima) = EMD.extrminmax(in)
@test minima == [3, 5, 7]
@test maxima == [2, 4, 6]
indzer = EMD.extrzeros(in)
@test indzer == [1, 2, 4, 5, 6, 7]

# now check the boundary conditions with interpolation, (hardcoding nsym)
(tmin, tmax, zmin, zmax) = EMD.boundarycheck(minima, maxima, t, in, in, 2)
@test tmin == [-1, 1, 3, 5, 7, 13, 15]
@test tmax == [-2, 0, 2, 4, 6, 10, 14]
@test zmin == [-2, -4, -2, -4, -3, -3, -4]
@test zmax == [0, 3, 3, 0, 1, 3, 1]

# testing mean/amp, which first requires boundary conditions, followed by
# extrema extraction.
(μ, nextr, nzer, amp) = EMD.meanamplitude(in)
@test isapprox(μ, [-0.229287791, 0.149086738, -0.312136628, -1.447260215, -1.959665698, -1.422545880, -0.568859012, 0.064862925, 0.410450924, 0.494524758], rtol=1e-8)
@test nextr == 6
@test nzer == 6
@test isapprox(amp, [3.770712209, 2.850913262, 1.687863372, 1.447260215, 2.040334302, 2.422545880, 2.431140988, 2.420602191, 2.455464774, 2.505475242], rtol=1e-8)



