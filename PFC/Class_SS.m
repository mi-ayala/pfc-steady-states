%Holds steady states properties to pass them around easily
classdef Class_SS < handle
	properties
		A; psi;
		E; E_error; E_intval;
		r_min; r_max;
		p_eig; n_eig; z_eig;
		unstable;
		eig_ok;
	end
end
