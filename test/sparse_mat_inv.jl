using Test
using IncompleteLU
using LinearAlgebra
using Hutchinson
using SparseArrays

@testset "Inverse of sparse matrices" begin
    
    A = sprand(10, 10, 0.2) + 10I
    @test abs(norm(sparse_mat_inv(A)) - norm(inv(Array(A)))) <= 1e-4

end