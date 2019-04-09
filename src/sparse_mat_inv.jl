
export sparse_mat_inv

using LinearAlgebra

#Breaking up U into DU
function gen_diag_mat(U::Matrix)

	l = length(U[:,1])
	D = zeros(l,l)
	for i = 1:l
		D[i,i] = U[i,i]
		U[i,i] = 1
		for j = i+1:l
			U[i,j] = U[i,j]/D[i,i]
		end
	end

	return D, U 
end

#Computing inverse of diagonal matrix
function diag_inv(D::Matrix)

	for i = 1:length(D[:,1])
		D[i,i] = 1/D[i,i]
	end

	return D
end

function sum_prod(A::Matrix, B::Matrix, i, j)
	sum = 0
	if i <= j
		for k = i+1:length(A[:,1])
			sum = sum + A[i,k]*B[k,j]
		end
	else
		for k = j+1:length(A[:,1])
			sum += A[i,k]*B[k,j]
		end
	end

	return sum
end

#Main Function
function sparse_mat_inv(A::Matrix)
	F = lu(A)
	L = F.L
	D, U = gen_diag_mat(F.U)
	D_inv = diag_inv(D)
	Z = Matrix{Float64}(I, length(A[:,1]), length(A[:,1]))
	for i = length(A[:,1]):-1:1

		#Diagonal
		Z[i,i] = D_inv[i,i] - sum_prod(U,Z,i,i)

		#Upper triangle
		for j = length(A[1,:]):-1:i+1
			if A[i,j] != 0
				Z[i,j] = - sum_prod(U,Z,i,j)
			end
		end

		#Lower triangle
		for j = i-1:-1:1
			if A[i,j] != 0
				Z[i,j] = - sum_prod(Z,L,i,j)
			end
		end

	end
	return Z

end

