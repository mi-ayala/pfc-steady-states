%Prepare the phase diagram from run data
compare_SS_folder = 'SS_PD/';

load(sprintf('%sSS_PD.mat', compare_SS_folder), 'req_pbar', 'req_beta');

%Initialize phase diagram data
pd_flag_bad_point = false(length(req_beta), length(req_pbar));
pd_by_id = NaN(size(pd_flag_bad_point)); pd_energy = pd_by_id;
pd_r_min = pd_by_id; pd_r_max = pd_by_id;
energy_constant = pd_by_id; energy_stripes = pd_by_id; energy_atoms = pd_by_id;

%Loop over the phase diagram
for cur_pbar_ind = 1:length(req_pbar)
	for cur_beta_ind = 1:length(req_beta)
		compare_SS_file = sprintf('%sSS_m%d_b%d/SS.mat', compare_SS_folder, cur_pbar_ind, cur_beta_ind);
		
		%Skip unavailable folders and leave them as NaN
		if ~exist(compare_SS_file, 'file')
			pd_by_id(cur_beta_ind, cur_pbar_ind) = 1;
			continue;
		end
		
		%Load the SS file
		load(compare_SS_file);
		
		%Determine the minimizer within the currently update candidates
		min_energy = list_SS{1}.E_intval;
		min_id = 1;
		for nn = 2:length(list_SS)
			%Ignore states with failed proofs
			if isnan(list_SS{nn}.r_min)
				pd_flag_bad_point(cur_beta_ind, cur_pbar_ind) = true;
				pd_flag_bad_point = true;
				continue;
			end
			
			%If the energy is bounded from the others, check if it is smaller
			cur_E = list_SS{nn}.E_intval;
			if emptyintersect(min_energy, cur_E)
				if cur_E < min_energy
					min_energy = cur_E;
					min_id = nn;
				end
			else
				%Attempt to remove translational symmetries for stripes and atoms
				if min_id == 2
					if abs(list_SS{min_id}.A(11,1) + list_SS{nn}.A(11,1)) > 10^-14
						fprintf('Symmetry failed for stripes.');
						pd_flag_bad_point(cur_beta_ind, cur_pbar_ind) = true;
					end
				elseif min_id == 3
					if abs(list_SS{min_id}.A(6,9) + list_SS{nn}.A(6,9)) > 10^-14
						fprintf('Symmetry failed for atoms.');
						pd_flag_bad_point(cur_beta_ind, cur_pbar_ind) = true;
					end
				else
					pd_flag_bad_point(cur_beta_ind, cur_pbar_ind) = true;
				end
			end
		end
		
		%Get some info about the ansatz
		energy_constant(cur_beta_ind, cur_pbar_ind) = list_SS{1}.E;
		energy_stripes(cur_beta_ind, cur_pbar_ind) = list_SS{2}.E;
		energy_atoms(cur_beta_ind, cur_pbar_ind) = list_SS{3}.E;
		
		%Get overall information
		pd_by_id(cur_beta_ind, cur_pbar_ind) = min_id;
		pd_energy(cur_beta_ind, cur_pbar_ind) = mid(min_energy);
		pd_r_min(cur_beta_ind, cur_pbar_ind) = list_SS{min_id}.r_min;
		pd_r_max(cur_beta_ind, cur_pbar_ind) = list_SS{min_id}.r_max;
	end
end

%Save the data and display the phase diagram
save(sprintf('%sSS_PD.mat', compare_SS_folder), 'req_pbar', 'req_beta', 'pd_by_id', ...
	'pd_flag_bad_point', 'pd_energy', 'pd_r_min', 'pd_r_max', 'energy_constant', 'energy_stripes', ...
	'energy_atoms');
PhaseDiagram_Show(req_pbar, req_beta, pd_by_id, 'RN_PhaseDiagram_Small.png')
