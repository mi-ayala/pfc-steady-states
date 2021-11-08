%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------MONTE CARLO STEADY STATES---------------------------------------
%%%%%%%%%%%%%% Test initial conditions randomly to find steady states with proofs
%%%%%%%%%%%%%% Gabriel Martine
%%%%%%%%%%%%%% Last updated: July 3 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; StartINTLAB();

%PFC and proof parameters
intval_psibar = intval('0.161060');
intval_beta = intval('0.07');
intval_nu = intval('1.01');

%Wavevector for two-mode PFC - use NaN for usual PFC
two_q = (2*cos(pi/5))^(-1);


%PFC and proof parameters
%intval_psibar = intval('0.5');
%intval_beta = intval('0.6');
%intval_nu = intval('1.01');

%Wavevector for two-mode PFC - use NaN for usual PFC
%two_q = NaN;

%Numerical parameters
N_atoms_x = 7;              %Number of atoms in the x direction
N_yx_ratio = 1;				%Controls the factor Ny = r*Nx, used for creating elongated domains
ppa = 64;					%Number of pixels per atoms
M = 60;						%Number of modes
compare_scale = 0.05;		%Initial noise scale

%Steady states file - records known coefficients and results into the following file
overwrite_SS = false;
compare_trials = 4;
compare_SS_file = 'SS/SS.mat';

%explicit_output_test = 'A.mat';		%Only used to test the output of Connect_SS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Array setup
A = zeros(M+1);
pfc_g = PFCGeneralData(two_q, M, N_atoms_x, N_yx_ratio, ppa, intval_psibar, intval_beta, intval_nu);

%When overwriting, remove the previous SS results if they exist
if overwrite_SS
	[~,~] = rmdir(fileparts(compare_SS_file), 's');
end

%Initialize the SS file if it does not exist yet
[~,~] = mkdir(fileparts(compare_SS_file));
if ~exist(compare_SS_file, 'file')
	list_SS = {}; list_counts = [];
	save(compare_SS_file, 'list_SS', 'list_counts', 'pfc_g');
else
	%Ensure the parameters have not changed since the last run
	old_pfc = getfield(load(compare_SS_file, 'pfc_g'), 'pfc_g');	
	if ~old_pfc.same_as(pfc_g)
		error('Numerical parameters changed, start a new run or overwrite!');
	end
	load(compare_SS_file);
end

%If testing a given output file, only use that initial condition
%if ~strcmp(explicit_output_test, '')
%	load(explicit_output_test);
%	[duplicate_id, list_SS, list_counts] = Accumulate_SS(output_A,pfc_g,list_SS,list_counts,compare_SS_file,false,false);
%	if isnan(duplicate_id)
%		error('Newton solver failed');
%	end
%	fprintf('\n');
%else
	%Otherwise, proceed through the trials and always test the four main ansatz at the end
	Atest = PrepareAnsatz(pfc_g, false);
	for nn = 1:compare_trials
		if nn > compare_trials-4
			A = Atest{compare_trials-nn+1};
		else
			A = compare_scale * randn(size(A));
			A(1,1) = pfc_g.psibar;
		end

		[duplicate_id, list_SS, list_counts] = Accumulate_SS(A,pfc_g,list_SS,list_counts,compare_SS_file,false,false);
		if isnan(duplicate_id)
			error('Newton solver failed');
		end
		fprintf('\n');
	end
%end

%Store the list file
save(compare_SS_file, '-append', 'list_SS', 'list_counts');

%Sort the steady states according to their energy
Sort_SS(compare_SS_file);
