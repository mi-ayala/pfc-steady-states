%Returns the coefficients of the four main ansatz
function ansatz = PrepareAnsatz(pfc_g, SKIP_DONUTS)	
	%Constant
	Ac = zeros(1+pfc_g.M);
	Ac(1,1) = pfc_g.psibar;
	
	%Stripes
	As = Ac;
	opt_amp = 2/sqrt(3)*sqrt(pfc_g.beta-3*pfc_g.psibar^2);
	if ~isreal(opt_amp)
		opt_amp = 0.0;
	end
	As(1+2*pfc_g.Ny, 1) = opt_amp/2;	%Using the + sign because it's the one that's stable!
	
	%Atoms
	Aa = Ac;
	opt_amp = -2/5*pfc_g.psibar - 2/sqrt(15)*sqrt(pfc_g.beta-12/5*pfc_g.psibar^2);
	if ~isreal(opt_amp)
		opt_amp = 0.0;
	end
	Aa(1+2*pfc_g.Ny, 1) = opt_amp/2; Aa(1+pfc_g.Ny, 1+pfc_g.Nx) = opt_amp/2;
	
	if SKIP_DONUTS
		%Pack everything
		ansatz = {Ac, As, Aa};
	else
		%Donuts
		Ad = Ac;
		opt_amp = -2/5*pfc_g.psibar + 2/sqrt(15)*sqrt(pfc_g.beta-12/5*pfc_g.psibar^2);
		if ~isreal(opt_amp)
			opt_amp = 0.0;
		end
		Ad(1+2*pfc_g.Ny, 1) = opt_amp/2; Ad(1+pfc_g.Ny, 1+pfc_g.Nx) = opt_amp/2;

		%Pack everything
		ansatz = {Ac, As, Aa, Ad};
	end
end
