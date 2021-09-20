%Creates a PFC grid that is oriented along the x-axis.
function pfc_g = PFCGeneralData(set_two_q, set_M, set_Nx, set_N_yx_ratio, set_ppa, ...
								set_intval_psibar, set_intval_beta, set_intval_nu)
	pfc_g = Class_PFCGeneralProperties;
	
	%Numerical grid data
	pfc_g.M = set_M;
	pfc_g.Nx = set_Nx;
	pfc_g.ppa = set_ppa;
	
	%Grid adapted to each lattice type
	pfc_g.two_q = set_two_q;
	if isnan(pfc_g.two_q)
		pfc_g.Ny = round(pfc_g.Nx/sqrt(3)) * set_N_yx_ratio;
		pfc_g.Lx = 4*pi/sqrt(3) * pfc_g.Nx; pfc_g.Ly = 4*pi * pfc_g.Ny;
		pfc_g.total_atoms = 2 * pfc_g.Nx * pfc_g.Ny;
	else
		pfc_g.Ny = round(pfc_g.Nx * set_N_yx_ratio);
		pfc_g.Lx = 2*pi * sqrt(2) * pfc_g.Nx; pfc_g.Ly = 2*pi * sqrt(2) * pfc_g.Ny;
		pfc_g.total_atoms = pfc_g.Nx * pfc_g.Ny;
	end
	
	%Representation grid
	pfc_g.x = linspace(0, pfc_g.Lx, pfc_g.Nx * pfc_g.ppa + 1); pfc_g.x(end) = [];
	pfc_g.y = linspace(0, pfc_g.Ly, pfc_g.Ny * pfc_g.ppa + 1); pfc_g.y(end) = [];
	
	%PFC parameters
	pfc_g.intval_psibar = set_intval_psibar; pfc_g.psibar = mid(pfc_g.intval_psibar);
	pfc_g.intval_beta = set_intval_beta; pfc_g.beta = mid(pfc_g.intval_beta);
	pfc_g.intval_nu = set_intval_nu; pfc_g.nu = mid(pfc_g.intval_nu);
	
	%Fourier matrices
	[pfc_g.lapl, pfc_g.gamma, pfc_g.nu_mat] = FourierMatrices(pfc_g.M, pfc_g, false);
	[pfc_g.intval_lapl, pfc_g.intval_gamma, pfc_g.intval_nu_mat] = FourierMatrices(pfc_g.M, pfc_g, true);
	
	%Constant energy
	pfc_g.constant_energy = ConstantEnergy(pfc_g.two_q, pfc_g.psibar, pfc_g.beta);
end
