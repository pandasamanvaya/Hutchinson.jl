module Hutchinson
using Hadamard
export inverse_mat_trace
include("sparse_mat_inv.jl")

#Generating random vectors
function random_vector(Z::AbstractMatrix)
	random = [-1 1]
	vec = zeros(length(Z[:,1]), 10)

	for i = 1 : length(Z[:,1])
		for j = 1 : 10
			vec[i,j] = random[rand(1:2)]
		end
	end

	return vec
end

#Main Function
function inverse_mat_trace(A::AbstractMatrix)

	Z = sparse_mat_inv(A)

	if length(Z[:,1])%4 == 0 || length(Z[:,1]) == 2
		H = hadamard(length(Z[:,1])) #using hadamard vectors
	else
		H = random_vector(Z) #using random vectors
	end
	
	trace = 0
	
	for i = 1 : length(H[1,:])
		e = H[:,i]
		trace += e'*Z*e
	end

	return (trace/length(H[1,:]))
end 

end # module
