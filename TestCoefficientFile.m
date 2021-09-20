%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------TEST COEFFICIENT FILE---------------------------------------
%%%%%%%%%%%%%% Tests a coefficient file in the appropriate PFC regime
%%%%%%%%%%%%%% Gabriel Martine
%%%%%%%%%%%%%% Last updated: Sep 18 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; StartINTLAB();

filename = 'CoefficientFiles/table_5/state_row1.mat';
load(filename)

% Radii polynomial proof results
[r_min, r_max, G] = RadiiPolyProof(A, pfc_g);

%Validated stability results
%[num_pos_eig, num_neg_eig, num_zero_eig, unstable_directions, _] = GetStability(G, true);

%Validated energy results
%[energy, energy_error] = RigorousEnergy(A, r_min, pfc_g);

%Show the PFC setup
fprintf('PFC proof:\n');
fprintf('intval (psibar, beta, nu): (%.2f, %.2f, %.2f)\n', pfc_g.intval_psibar.mid, pfc_g.intval_beta.mid, pfc_g.intval_nu.mid);
if ~isnan(pfc_g.two_q)
	fprintf('Two-mode PFC with q: %.3f\n', pfc_g.two_q);
end
fprintf('domain (N_x, N_y): (%.2f, %.2f)\n', pfc_g.Nx, pfc_g.Ny);
fprintf('ppa: %.2f\n', pfc_g.ppa);
fprintf('M: %d\n', pfc_g.M);
