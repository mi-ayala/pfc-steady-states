%Derivative of F for 2D PFC with periodic Neumann boundary conditions
function mat_DF = DF_2D(A, pfc_g, USE_INTLAB)
	M = pfc_g.M;
	use_lapl = pfc_g.lapl;
	use_gamma = pfc_g.gamma;
	if USE_INTLAB
		if ~isa(A, 'intval') error('A is not an intval!'); end
		use_lapl = pfc_g.intval_lapl;
		use_gamma = pfc_g.intval_gamma;
	end
	
	%The derivative of the cubic term gives 3x(four/two/one) terms depending on sigma, tau
	AA = ConvAk_2D(A, 2);	
	
	%Computation of the derivative matrix (already in matrix form)
	mat_DF = zeros((1+M)^2);
	if USE_INTLAB
		mat_DF = intval(mat_DF);
	end
	
	%DF's column 0 corresponds to the partial derivatives d_0,0
	mat_DF(:, 1) = SwitchMat2Vec(AA(1+(0:M), 1+(0:M)));
	
	%DF's columns 1 through M+1 correspond to the partial derivatives d_sigma,0
	for cur_sigma = 1:M
		mat_DF(:, 1+cur_sigma) = SwitchMat2Vec(AA(1+(0:M), 1+abs(cur_sigma + (0:M))) + AA(1+(0:M), 1+abs(cur_sigma - (0:M))));
	end
	
	for cur_tau = 1:M
		fac_plus = 1+abs(cur_tau + (0:M));
		fac_minus = 1+abs(cur_tau - (0:M));
		
		%Other matrix columns, Column 0
		mat_DF(:, 1+(M+1)*cur_tau) = SwitchMat2Vec(AA(fac_plus, 1+(0:M)) + AA(fac_minus, 1+(0:M)));
		
		%General case
		for cur_sigma = 1:M
			mat_DF(:, 1+cur_sigma+(M+1)*cur_tau) = SwitchMat2Vec(AA(fac_plus, 1+abs(cur_sigma + (0:M))) + ...
				AA(fac_plus, 1+abs(cur_sigma - (0:M))) + AA(fac_minus, 1+abs(cur_sigma + (0:M))) + ...
				AA(fac_minus, 1+abs(cur_sigma - (0:M))));
		end
	end
	
	%Add the gamma linear terms
	mat_DF = 3*mat_DF + diag(SwitchMat2Vec(use_gamma));
	
	%Multiply by the Laplacian coefficients and fix the j=k=0 coefficients (mat_DF(1, :) = 0; is already fixed)
	mat_DF = SwitchMat2Vec(use_lapl) .* mat_DF;
	mat_DF(1, 1) = 1;
end
