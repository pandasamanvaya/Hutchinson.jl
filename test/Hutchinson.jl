using Test
using SparseArrays
using LinearAlgebra
using Hutchinson

@testset "trace of matrix inverse" begin
    A = sprand(7, 7, 0.2) + 7I
    @test abs(inverse_mat_trace(A) - tr(inv(A))) <= 0.01
end