%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------TEST ANSATZ---------------------------------------
%%%%%%%%%%%%%% Test the four main ansatz (from simple coefficients)
%%%%%%%%%%%%%% Gabriel Martine
%%%%%%%%%%%%%% Last updated: Sep 19 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;

% Intval must be available through INTLAB - it can be obtained at https://www.tuhh.de/ti3/rump/intlab/
% StartINTLAB();

% The current project must be loaded in the path
% addpath(genpath('.')

params_psibar = [0.07, 0.3, 0.5]; params_beta = [0.025, 0.5, 1.0]; params_M = [20, 30, 40];

for mm = 1:length(params_psibar)
	pfc_g = PFCGeneralData(NaN, params_M(mm), 4, 1, 128, intval(params_psibar(mm)), intval(params_beta(mm)), intval(1.05));
	Atest = PrepareAnsatz(pfc_g, false);
	
	%Test the ansatz
	for nn = 1:4
		switch nn
			case 1
				fprintf('Testing Constant state\n')
			case 2
				fprintf('Testing Stripes state\n')
			case 3
				fprintf('Testing Atoms state\n')
			case 4
				fprintf('Testing Donuts state\n')
		end
		Afinal = NewtonSolver_2D(Atest{nn}, pfc_g, 100, 10^-16, false);
		[r_min, r_max, G_PFC] = RadiiPolyProof(Afinal, pfc_g);
		[E, Ed] = RigorousEnergy(Afinal, r_min, pfc_g);
		[p_eig, ~, ~, ~, ~] = GetStability(G_PFC, false);
		
		%Output info
		ss_str = sprintf('&$(%.3f, %.3f)$&\\makecell{\\vspace*{0.05cm}$%.2f$\\\\$%.1e}$}&\\makecell{\\vspace*{0.05cm}$%.1e}$\\\\$%.1e}$}&\\makecell{$%.3e}$\\\\$%.1e}$}&%d\\\\', ...
			pfc_g.psibar, pfc_g.beta, NormMatrixNu(Afinal, pfc_g.nu_mat), NormMatrixNu(Afinal-Atest{nn}, pfc_g.nu_mat), r_min, r_max, E-pfc_g.constant_energy, Ed, p_eig-1);
		ss_str = strrep(ss_str, '0.0e+00}', '<\epsilon');
		ss_str = strrep(ss_str, '0.000e+00}', '<\epsilon');
		ss_str = strrep(ss_str, 'e-0', '\cdot10^{-');
		ss_str = strrep(ss_str, 'e-', '\cdot10^{-');
		fprintf('\\hline\n%s\n\n\n', ss_str);
		
		%Visualize
		%SavePhaseImage(GetPhase(Afinal, pfc_g), pfc_g, sprintf('Im%d.png', nn));
		%pause(1.0);
	end
end
