%Holds grid properties to pass them around easily
classdef Class_PFCGeneralProperties < handle
	properties
		%Numerical grid properties
		M;
		Nx; Ny; Lx; Ly;
		ppa; total_atoms;
		x; y;
		
		%PFC parameters
		two_q;
		psibar; beta; nu;
		intval_psibar; intval_beta; intval_nu;
		
		%Fourier matrices
		lapl; gamma; nu_mat;
		intval_lapl; intval_gamma; intval_nu_mat;
		
		%Constant energy
		constant_energy;
	end
	
	methods
		%Ensures that this object contains the same properties as another
		function v = same_as(self, other)
			assert(isa(other, 'Class_PFCGeneralProperties'));
			v = false;
			if (isnan(self.two_q) && isnan(other.two_q)) || (self.two_q - other.two_q) < 10^-10
				v = isequaln([self.M, self.Nx, self.Ny, self.ppa, self.psibar, self.beta, self.nu], ...
					[other.M, other.Nx, other.Ny, other.ppa, other.psibar, other.beta, other.nu]);
			end
		end
	end
end