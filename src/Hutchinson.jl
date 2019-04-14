module Hutchinson
using Hadamard
export inverse_mat_trace
include("sparse_mat_inv.jl")

function hadamard_vector(Z::AbstractMatrix)
	H = hadamard(length(Z[:,1]))
	trace = 0
	for i = 1 : length(H[1,:])
		e = H[:,i]
		trace += e'*Z*e
	end
	return (trace/length(H[1,:]))
end

function random_vector(Z::AbstractMatrix)
	random = [-1 1]
	e = zeros(length(Z[:,1]))

	for i = 1 : length(Z[:,1])
		e[i] = random[rand(1:2)]
	end

	return e'*Z*e

function inverse_mat_trace(A::AbstractMatrix)

	Z = sparse_mat_inv(A)
	if length(Z[:,1])%4 == 0 || length(Z[:,1]) == 2
		trace = hadamard_vector(Z)
	else
		trace = random_vector(Z)
	end
	
	return trace
end 

end # module
