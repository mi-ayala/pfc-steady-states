%Function F for 2D PFC with periodic Neumann boundary conditions
function vec_F = F_2D(A, pfc_g, USE_INTLAB)
	use_psibar = pfc_g.psibar;
	use_lapl = pfc_g.lapl;
	use_gamma = pfc_g.gamma;
	if USE_INTLAB
		if ~isa(A, 'intval') error('A is not an intval!'); end
		use_psibar = pfc_g.intval_psibar;
		use_lapl = pfc_g.intval_lapl;
		use_gamma = pfc_g.intval_gamma;
	end
	
	%Most of F is given by laplacian_jk (gamma_jk A_jk + (A*A*A)_jk)
	AAA = ConvAk_2D(A, 3);
	mat_F = use_lapl .* (use_gamma .* A + AAA(1+(0:pfc_g.M), 1+(0:pfc_g.M)));
	
	%The first coefficient is a_00 - psibar with derivative delta_0 delta_0
	mat_F(1,1) = A(1,1) - use_psibar;
	
	%Reshape the matrix into its vector formulation
	vec_F = SwitchMat2Vec(mat_F);
end
