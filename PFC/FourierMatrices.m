%Computes the fourier matrices for the Laplacian, gamma and nu_mat
function [lapl, gamma, nu_mat] = FourierMatrices(M, pfc_g, USE_INTVAL)
	[JJ, KK] = meshgrid(0:M);
	
	%Compute the Fourier Laplacian
	use_beta = pfc_g.beta;
	use_nu = pfc_g.nu;
	if USE_INTVAL
		use_beta = pfc_g.intval_beta;
		use_nu = pfc_g.intval_nu;
		lapl = -((intval('2')*pi/pfc_g.Lx * JJ).^2 + (intval('2')*pi/pfc_g.Ly * KK).^2);
	else
		lapl = -((2*pi/pfc_g.Lx * JJ).^2 + (2*pi/pfc_g.Ly * KK).^2);
	end
	
	%Linear term of the PDE using PFC or two-mode PFC
	if isnan(pfc_g.two_q)
		gamma = (lapl + 1).^2 - use_beta;
	else
		gamma = (lapl + 1).^2 .* (lapl + pfc_g.two_q^2).^2 - use_beta;
	end
	
	%Coefficient matrix for l1nu computations
	nu_mat = WeightMatrix(M) .* use_nu .^ (JJ + KK);
end
