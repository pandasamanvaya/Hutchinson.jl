using Test
using SparseArrays
using LinearAlgebra
using Hutchinson

@testset "trace of matrix inverse" begin
    A = sprand(7, 7, 0.2) + 7I
    @test abs(inverse_mat_trace(A) - tr(inv(Array(A)))) <= 0.01

    A = sprand(12, 12, 0.4) + 12I
    @test abs(inverse_mat_trace(A) - tr(inv(Array(A)))) <= 0.01
end