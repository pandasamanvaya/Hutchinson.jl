using Hutchinson
using LinearAlgebra
using Test

A = [4 0 0 0; 0 1 1 0; 0 0 2 0; 3 0 0 1]
A_inv = inv(A)
trace = tr(A_inv)

@test inverse_mat_trace(A) == trace 