load("extras/mincg.jl")

function myf1(X::Matrix)
       return (sin(X[1]) + cos(X[2]) , [cos(X[1]) -sin(X[2])]')
 end

(minX, steps, i) =  mincg(myf1, [1 1]', 100)
(minf1, minD1 ) = myf1(minX)

println("Minimum of $minf1 @ $(minX[1]), $(minX[2]) in $i Iterations")

@assert_approx_eq minX[1]  -pi/2
@assert_approx_eq minX[2]  pi
@assert abs(-2 - minf1)  < 2*realmin()


function myf2(X::Matrix)
	c=cos(X[1]*X[2])
	return (sin(X[1]*X[2]), [ X[2]*c  X[1]*c ]')
end

(minX, steps, i) =  mincg(myf2, [-1 1]', 100)
(minf1, minD1 ) = myf2(minX)
println("Minimum of $minf1 @ $(minX[1]), $(minX[2]) in $i Iterations")

# @assert_approx_eq minX[1]  -1.37517
# @assert_approx_eq minX[2]  1.14226
@assert abs(-1 - minf1)  < 2*realmin()