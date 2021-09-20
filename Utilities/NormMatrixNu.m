%Computes the l1nu(N^2) norm of a coefficient matrix
function nnu = NormMatrixNu(A, nu_matrix)
	nnu = sum(sum(abs(A) .* nu_matrix));
end
