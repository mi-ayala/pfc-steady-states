%Computes the B(l1nu(N^2)) norm of a 4-tensor in its matrix formulation
function nnu = NormTensorNu(Q, nu_matrix)
	nu_vector = SwitchMat2Vec(nu_matrix);
	nnu = max((abs(Q)*nu_vector) ./ nu_vector);
end
