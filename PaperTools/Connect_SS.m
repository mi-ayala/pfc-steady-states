%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------CONNECT STEADY STATES-------------------------------------------
%%%%%%%%%%%%%% Test initial conditions randomly to capture steady states and prove them.
%%%%%%%%%%%%%% 1D unstable manifold are simple, 2D can be managed using a proper discretization of the circle.
%%%%%%%%%%%%%% Higher order steady states are too complicated for now. Attach pictures instead of labels when rendering.
%%%%%%%%%%%%%% Gabriel Martine
%%%%%%%%%%%%%% Last updated: Nov 18 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; StartINTLAB();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PARAMETERS
compare_SS_file = 'SS/SS.mat';
output_A_file = 'A.mat';
state_id = 9;
override_test_directions = [1];		%Leave blank for automatic testing, otherwise specify the unstable directions
render_thesis = false;
auto_montecarlo = true;				%Set to true to test invalid flow endpoints

%Perturbation
test_radius = 0.001;
max_iterations = 2600;
check_flow_step = 1;			%How often to check convergence of the flow. Be careful, too large and we might miss states
energy_steps = 4;

%Stop when given tolerances are reached
energy_tolerance = 1e-5;
phase_tolerance = 3e-2;			%WATCHA CHANGE THIS TO 1E-2

%2D unstable manifold parameters
N_tests = 100;					%Sweep is CCW
focus_angle = atan2(0, 1);		%atan2(Y, X), Y=0, X=1
focus_spread = 1.0*(2*pi);

%PFC flow parameters
tau = 0.5;
C = 0.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load the SS file
load(compare_SS_file);
auto_montecarlo_new_states = false;

%Check if overriding test directions
if isempty(override_test_directions)
	%Stop if the input is bad or there are zero eigenvalues
	if list_SS{state_id}.z_eig > 0
		fprintf('Unable to continue, zero eigenvalue(s) detected for this steady state!\n');
		return;
	elseif list_SS{state_id}.p_eig == 1
		fprintf('Unable to continue, the steady state is stable!\n');
		return;
	elseif list_SS{state_id}.p_eig > 3
		fprintf('Higher order unstable manifolds are not implemented yet!\n');
		return;
	end
	
	%First we must remove the trivial direction associated to psibar (its (1,1) coefficient is MOST LIKELY the greatest!)
	unstable_dirs = list_SS{state_id}.unstable;
	[~, remove_dir] = max(abs(cellfun(@(obj)obj(1,1), unstable_dirs)));
	unstable_dirs(remove_dir) = [];
else
	unstable_dirs = {};
	for nn = 1:min(2, length(override_test_directions))
		unstable_dirs{end+1} = list_SS{state_id}.unstable{override_test_directions(nn)};
	end
end

%--------------TEST PHASES-------------------------------------------
%Compute the real space unstable directions in the span of cosines
psi_unstable_dir = {};
for cu = 1:length(unstable_dirs)
	psi_unstable_dir{end+1} = GetPhase(unstable_dirs{cu}, pfc_g);
end

%FOR SPECIAL CASE STARTING AT SMALL HEX -> GOES TO UNSTABLE HORIZONTAL INTERMEDIATE
%psi_unstable_dir{1} = GetPhase(0.68*list_SS{state_id}.unstable{4}+0.24*list_SS{state_id}.unstable{7}-0.08*list_SS{state_id}.unstable{5}, pfc_g);

%FOR SPECIAL CASE STARTING AT BIG HEX -> GOES TO STABLE HEXAGONAL PATCH
%psi_unstable_dir{1} = GetPhase(-0.23*list_SS{state_id}.unstable{1}+0.26*list_SS{state_id}.unstable{13}+0.26*list_SS{state_id}.unstable{15}-0.25*list_SS{state_id}.unstable{16}, pfc_g);

%Prepare perturbations. Higher dimensions need to discretize the unit sphere about a focus angle
switch length(psi_unstable_dir)
	case 1	%Should be [1, -1] for correct plotting
		test_vectors = [1]; fprintf('WARNING, DISABLING BOTH DIRECTIONS.\n');
	case 2
		angles = focus_angle + focus_spread/N_tests * (-floor(N_tests/2):floor(N_tests/2));
		N_tests = length(angles);
		test_vectors = [cos(angles)', sin(angles)'];
	otherwise
		return;
end


%--------------FLOW PREPARATION-------------------------------------------
%PFC grid and Fourier matrices
hx = pfc_g.x(2)-pfc_g.x(1); hy = pfc_g.y(2)-pfc_g.y(1); grid_N = pfc_g.ppa * [pfc_g.Nx, pfc_g.Ny];
[X, Y] = meshgrid(0:(grid_N(1)-1), 0:(grid_N(2)-1));
lapl_PFC = 2.0/hx^2 * (cos(2.0*pi*X/grid_N(1)) - 1.0) + 2.0/hy^2 * (cos(2.0*pi*Y/grid_N(2)) - 1.0);
G_PFC = (1.0 - tau * lapl_PFC .* ((lapl_PFC+1.0).^2 + C - pfc_g.beta)).^(-1);


%--------------ITERATION OVER ALL TEST DIRECTIONS-------------------------------------------
%Iterate over all test directions
connections = zeros(length(test_vectors), 1); energy_flows = {};
for tv = 1:length(test_vectors)
	fprintf('Flowing test direction %d.\n', tv);
	
	%Compute the starting point for the flow using the current test direction
	switch length(psi_unstable_dir)
		case 1
			psi = test_vectors(tv) * psi_unstable_dir{1};
		case 2
			psi = test_vectors(tv, 1) * psi_unstable_dir{1} + test_vectors(tv, 2) * psi_unstable_dir{2};
	end
	psi = list_SS{state_id}.psi + test_radius * psi;
	fft_psi = fft2(psi);
	
	%Determine the steady states allowed for comparison by disregarding those with higher or equal energy
	cur_energy_flow = list_SS{state_id}.E;
	compare_list = 1:length(list_SS);
	compare_list(cellfun(@(obj)obj.E >= cur_energy_flow(1), list_SS)) = [];
	
	%--------------PFC FLOW-------------------------------------------
	for iter = 1:max_iterations
		%Compute the current PFC energy (in terms of the phase)
		cur_energy_flow(end+1) = mean2(0.5*psi .* real(ifft2((lapl_PFC+1.0).^2 .* fft2(psi))) + ...
			0.25*(psi.^2 - pfc_g.beta).^2);
				
		%Check if the energy is close enough to that of other steady states
		if ~mod(iter, check_flow_step)
			for compare_id = compare_list
				if abs(cur_energy_flow(end) - list_SS{compare_id}.E) < energy_tolerance
					%Use a normalized L1 error
					l1_error = mean2(abs(psi - list_SS{compare_id}.psi));
					if l1_error < phase_tolerance
						fprintf('Flow near state %d. Stopping after %d iterations.\n', compare_id, iter);
						connections(tv) = compare_id;
						break;
					end
				end
			end
			
			if connections(tv) ~= 0
				break;
			end
		end
		
		%PFC Flow
		fft_psi = G_PFC .* (fft_psi + tau * lapl_PFC .* (fft2(psi.^3) - C*fft_psi));
		psi = real(ifft2(fft_psi));
		fft_psi = fft2(psi);
		
		if ismember(iter, 4*[341, 416, 486, 563])
			SavePhaseImage(psi, pfc_g, sprintf('img_%d.png', iter));
		end
		
		
		
		%Only plot every few steps
		if ~mod(iter, 8)
			figure(2); imagesc(psi); colormap(flipud(bone));
			title(sprintf('Step %d, Average (abs) offset phase: %.1e', iter, mean2(abs(psi-pfc_g.psibar))));
			pause(0.01);
		end
	end
	energy_flows{end+1} = cur_energy_flow(1:energy_steps:length(cur_energy_flow));
	
	%--------------OUTPUT FOR THE CURRENT DIRECTION-------------------------------------------
	%If no connecting state has been found, print a sad message
	if connections(tv) == 0
		fprintf('Flow did not get close to another steady state. ');
		output_A = GetCoefficients(psi, pfc_g.M, pfc_g, false);
		
		%Verify the new state if needed
		if auto_montecarlo
			fprintf('Attempting to verify the flow endpoint.\n');
			load(compare_SS_file);
			[duplicate_id, list_SS, list_counts] = Accumulate_SS(...
				output_A,pfc_g,list_SS,list_counts,compare_SS_file,false,false);
			save(compare_SS_file, '-append', 'list_SS', 'list_counts');
			auto_montecarlo_new_states = true;
			connections(tv) = duplicate_id;
		else
			fprintf('Saving the flow endpoint.\n');
			save(output_A_file, 'output_A');
		end
	end
end


%--------------VISUALIZATION-------------------------------------------
%REVIEW THIS SECTION JUST TO MAKE SURE!!!

switch length(psi_unstable_dir)
	case 1
		figure(3); hold on; xlabel('Iteration'); ylabel('Energy');
		
		%Plot the starting energy
		scatter(0, list_SS{state_id}.E, 25, 'k', 'filled');
		text(0, list_SS{state_id}.E, {sprintf('%d', state_id), '', ''});
		
		%Forward direction
		plot(1:length(energy_flows{1}), energy_flows{1}, 'b', 'linewidth', 2.0);
		text(length(energy_flows{1}), energy_flows{1}(end), sprintf('%d', connections(1)));
		
		%Backward direction
		plot(-1:-1:-length(energy_flows{2}), energy_flows{2}, 'b', 'linewidth', 2.0);
		text(-length(energy_flows{2}), energy_flows{2}(end), sprintf('%d', connections(2)), ...
			'horizontalalignment', 'right');
	
	case 2
		figure(3); hold on; xlabel('first direction $(0, 2N_y)$', 'Interpreter', 'latex');
		ylabel('Energy', 'Interpreter', 'latex'); zlabel('Energy');
		colormap('jet');
		
		%Create line group so they can all be removed
		g_line = hggroup();
		
		%Plot each direction as a straight line going away from the origin
		for cang = 1:N_tests
			xx = cos(angles(cang))*(1:length(energy_flows{cang}));
			yy = sin(angles(cang))*(1:length(energy_flows{cang}));
			scatter3(xx, yy, energy_flows{cang}, 25, energy_flows{cang}, 'filled');
			
			if ~render_thesis
				plot3(xx, yy, energy_flows{cang}, 'k', 'linewidth', 1.0, 'parent', g_line);
				
				%Reducing text usage by omitting repeated states
				if (cang == 1) || (cang == N_tests) || ...
						(connections(cang) ~= connections(cang-1)) || (connections(cang) ~= connections(cang+1))
					text(xx(1)+xx(end), yy(1)+yy(end), energy_flows{cang}(end), sprintf('%d', connections(cang)), ...
						'FontWeight', 'bold');
				end
			end
		end
		
		%Plot the starting energy
		scatter3(0, 0, list_SS{state_id}.E, 25, 'k', 'filled');
		if ~render_thesis
			text(0, 0, 0, sprintf('%d', state_id));
		end
end

if render_thesis
	%Make the graph a bit nicer
	xticklabels([]); yticklabels([]); xticks([]); yticks([]); box on; axis equal;
	set(findall(gcf, '-property', 'FontSize'), 'FontSize', 16)
	ax = gca; outerpos = ax.OuterPosition; ti = ax.TightInset;
	ax.Position = [outerpos(1)+ti(1), outerpos(2)+ti(2)+0.01, outerpos(3)-ti(1)-ti(3), outerpos(4)-ti(2)-ti(4)-0.01];
	saveas(2, 'out.png');
end

%Sort the new steady states
if auto_montecarlo_new_states
	Sort_SS(compare_SS_file);
end

%--------------OUTPUT-------------------------------------------
%Now that all test directions have been explored, output the result
connections = unique(connections);
if ismember(0, connections)
	fprintf('Some test directions failed to converge!\n');
	connections(connections == 0) = [];
end
fprintf('Tests complete. State %d connects to steady state(s): %s.\n', state_id, sprintf('%d ', connections));
