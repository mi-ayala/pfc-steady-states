%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------PHASE DIAGRAM---------------------------------------
%%%%%%%%%%%%%% Prepares a phase diagram in (psibar, beta) by looping through the inputs
%%%%%%%%%%%%%% Gabriel Martine
%%%%%%%%%%%%%% Last updated: July 7 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; StartINTLAB();

%This is for localized regime
%PFC Parameters, ignoring the beta = 0 case.
req_pbar = linspace(0.4, 0.6, 28);
req_beta = linspace(0.5, 0.7, 28);
multiscale_inclusion = @(m,b) (b>0.1+0.72*37/15*(m)^2) && (b<-0.18+1.5*37/15*(m)^2);

%Proof parameters
proof_nu = intval('1.01');
Nx = 4;
M = 30;


%THIS IS FOR THE SAMLL PHASE DIAGRAM
% %PFC Parameters, ignoring the beta = 0 case.
% req_pbar = linspace(0, 0.1, 10);
% req_beta = linspace(0, 0.05, 10);
% 
% %Proof parameters
% proof_nu = intval('1.01');
% Nx = 4;
% M = 20;

%Other parameters
ppa = 64;
compare_scale = 0.05;		%Initial noise scale
compare_trials = 14;
compare_SS_folder = 'SS_PD/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create the directory
[~,~] = mkdir(compare_SS_folder);

%Save requested parameters
save(sprintf('%sSS_PD.mat', compare_SS_folder), 'req_pbar', 'req_beta');

%Create a placeholder state, assuming the constant state is always obtained
A = zeros(M+1);
placeholder_SS = Class_SS; placeholder_SS.E_intval = intval(10^10); placeholder_SS.E = NaN;
placeholder_SS.A = NaN*A; placeholder_SS.r_max = NaN; placeholder_SS.p_eig = NaN;

%Loop over the phase diagram
for cur_beta_ind = 1:length(req_beta)
	for cur_pbar_ind = 1:length(req_pbar)
		%Use the global result, or multiscale inclusions, to do less work
		if (req_beta(cur_beta_ind) <= req_pbar(cur_pbar_ind)^2) || ...
				~multiscale_inclusion(req_pbar(cur_pbar_ind), req_beta(cur_beta_ind))
			continue;
		end
		
		%Initialize PFC data
		pfc_g = PFCGeneralData(NaN, M, Nx, 1, ppa, intval(num2str(req_pbar(cur_pbar_ind))), ...
			intval(num2str(req_beta(cur_beta_ind))), proof_nu);
				
		%Initialize the SS file if it does not exist yet
		compare_SS_file = sprintf('%sSS_m%d_b%d/SS.mat', compare_SS_folder, cur_pbar_ind, cur_beta_ind);
		[~,~] = mkdir(fileparts(compare_SS_file));
		trial_nn = 1:compare_trials;
		Atest = [];
		if ~exist(compare_SS_file, 'file')
			list_SS = {}; list_counts = [];
			Atest = PrepareAnsatz(pfc_g, true);
		else
			%Skip if already computed
			continue;
		end
		
		%Proceed through the trials by starting with the known ansatz (CONSTANT, STRIPES, ATOMS)
		fprintf('------------TRIAL (%d, %d): (psibar, beta) = (%.2f, %.2f)------------\n', ...
			cur_pbar_ind, cur_beta_ind, pfc_g.psibar, pfc_g.beta);
		for nn = trial_nn
			%Always test the four main ansatz at the start
			if nn < 4
				A = Atest{nn};
			else
				A = compare_scale * randn(size(A));
				A(1,1) = pfc_g.psibar;
			end
			
			%Stop looping if the Newton solver fails. Inneficient to save and reload all the time, but somewhat simpler
			[duplicate_id, list_SS, list_counts] = Accumulate_SS(A,pfc_g,list_SS,list_counts,compare_SS_file,false,false);
			if isnan(duplicate_id)
				fprintf('Newton solver failed');
			end
			
			%Check if the stripes or atoms state have failed and add the placeholder garbage instead
			if (nn == 2 || nn == 3) && (duplicate_id ~= nn)
				list_SS{end+1} = placeholder_SS; list_counts(end+1) = 0;
			end
			
			%fprintf('\n');
		end
		
		%Store the list file
		save(compare_SS_file, 'list_SS', 'list_counts', 'pfc_g');
	end
end
