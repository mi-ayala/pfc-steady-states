%I disabled stability detection for the phase diagram

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------ACCUMULATE STEADY STATES-----------------------------------------
%%%%%%%%%%%%%% Adds a new steady state of PFC (with periodic Neumann boundary conditions) to the comparison file
%%%%%%%%%%%%%% The steady state is verified rigorously using the radii polynomial approach and interval arithmetic
%%%%%%%%%%%%%% Unstable directions are found, and the information is saved to the comparison file
%%%%%%%%%%%%%% Gabriel Martine
%%%%%%%%%%%%%% Last updated: July 3 2019

%--------------INPUT------------------------------------------------------------
%A								%Starting guess
%pfc_g, list_SS, list_counts	%Main data set
%verify_eigenvalues				%Skip the verification of unstable directions
%PRINT_TIMING					%Whether to display the time taken in each operation

%--------------OUTPUT-----------------------------------------------------------
%The ID of the steady state, may be duplicated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [duplicate_id, list_SS, list_counts] = Accumulate_SS(A, pfc_g, list_SS, list_counts, compare_SS_file, ...
	verify_eigenvalues, PRINT_TIMING)

%--------------NEWTON SOLVER AND DUPLICATE REMOVAL------------------------------
tic;															%Use 10^-13 if the solver fails to converge
[A, failed_to_converge, iteration_step] = NewtonSolver_2D(A, pfc_g, 300, 10^-13, false);
if failed_to_converge
	duplicate_id = NaN;
	printf('WARNING!!! Newton solver failed to converge!\n');
	return;
else
	fprintf('Solver converged in %d steps.', iteration_step);
end

%Check if the current solution is within the proof radius of the known steady states
for nn = 1:length(list_SS)
	if NormMatrixNu(list_SS{nn}.A - A, pfc_g.nu_mat) <= list_SS{nn}.r_max
		duplicate_id = nn;
		list_counts(duplicate_id) = list_counts(duplicate_id) + 1;
		fprintf('Detected duplicate SS %d.\n', duplicate_id);
		return;
	end
end
if PRINT_TIMING fprintf('Time elapsed in the Newton solver: %.2fs\n', toc); end


%--------------PROOF, ENERGY, EIGENVALUE CHECK AND PLOT-------------------------
tic;
[r_min, r_max, G_PFC] = RadiiPolyProof(A, pfc_g);
if PRINT_TIMING fprintf('Time elapsed in the proof: %.2fs\n', toc); end

tic;
[energy, energy_error] = RigorousEnergy(A, r_min, pfc_g);
if PRINT_TIMING fprintf('Time elapsed computing the energy: %.2fs\n', toc); end

tic;
[num_pos_eig, num_neg_eig, num_zero_eig, unstable_directions, eigenvalues_ok] = GetStability(G_PFC, verify_eigenvalues);
if PRINT_TIMING fprintf('Time elapsed computing unstable directions: %.2fs\n', toc); end


%--------------ADDITION TO THE SS FILE (not necessarily verified!)---------
list_SS{end+1} = Class_SS; list_counts(end+1) = 1;
list_SS{end}.A = A; list_SS{end}.psi = GetPhase(A, pfc_g);
list_SS{end}.E = energy; list_SS{end}.E_error = energy_error; list_SS{end}.E_intval = midrad(energy, energy_error);
list_SS{end}.r_min = r_min; list_SS{end}.r_max = r_max;
list_SS{end}.p_eig = num_pos_eig; list_SS{end}.n_eig = num_neg_eig; list_SS{end}.z_eig = num_zero_eig;
list_SS{end}.unstable = unstable_directions; list_SS{end}.eig_ok = eigenvalues_ok;

%Duplicate SS is the new number of steady states
duplicate_id = length(list_SS);

%Plot the steady state
figure('visible','off');
SavePhaseImage(list_SS{duplicate_id}.psi, pfc_g, [fileparts(compare_SS_file), sprintf('/SS_%d.png', duplicate_id)]);
end
