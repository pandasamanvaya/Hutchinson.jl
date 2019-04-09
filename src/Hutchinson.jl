#module Hutchinson
using Hadamard
export inverse_mat_trace
include("sparse_mat_inv.jl")

function inverse_mat_trace(A::AbstractMatrix)
	Z = sparse_mat_inv(A)

	H = hadamard(length(A[:,1]))
	trace = 0
	for i = 1 : length(H[1,:])
		e = H[:,i]
		trace += e'*Z*e
	end
	return (trace/length(H[1,:]))
end 

A = [4 0 0 0; 0 1 1 0; 0 0 2 0; 3 0 0 1]
println(inverse_mat_trace(A))
#end # module
