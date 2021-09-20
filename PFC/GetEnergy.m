%Computes the PFC energy of a coefficient matrix
function sol_energy = GetEnergy(A, pfc_g, USE_INTLAB)
	use_lapl = pfc_g.lapl;
	beta_mat = 0*A;
	if USE_INTLAB
		if ~isa(A, 'intval') error('A is not an intval!'); end
		use_lapl = pfc_g.intval_lapl;
		beta_mat = intval(beta_mat);
		beta_mat(1,1) = pfc_g.intval_beta;
	else
		beta_mat(1,1) = pfc_g.beta;
	end
	
	%Wave term of the energy, for PFC or multimode PFC 
	if isnan(pfc_g.two_q)
		wave_matrix = ConvAk_2D((use_lapl + 1) .* A, 2);
	else
		wave_matrix = ConvAk_2D((use_lapl + 1) .* (use_lapl + pfc_g.two_q^2) .* A, 2);
	end
	
	%Well term of the energy
	well_matrix = ConvAk_2D(ConvAB_2D(A, A) - beta_mat, 2);
		
	%The integral of the Fourier series returns the first component only
	sol_energy = 0.5*wave_matrix(1,1) + 0.25*well_matrix(1,1);
end
