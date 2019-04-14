using Test
using IncompleteLU
using LinearAlgebra
using Hutchinson
using SparseArrays

@testset "Inverse of sparse matrices" begin
    
    A = sprand(10, 10, 0.2) + 10I
    @test abs(norm(sparse_mat_inv(A)) - norm(inv(Array(A)))) <= 1e-4

    A = sprand(12, 12, 0.8) + 10I
    @test_throws ArgumentError sparse_mat_inv(A)
end