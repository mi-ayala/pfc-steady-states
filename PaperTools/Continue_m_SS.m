%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------CONTINUE STEADY STATES------------------------------------------
%%%%%%%%%%%%%% This follows a branch using pseudo-arclength continuation. The branch starts at some value of
%%%%%%%%%%%%%% (psibar, beta) for an SS_file with the initial A matrix indicated by the given state_id
%%%%%%%%%%%%%% Gabriel Martine
%%%%%%%%%%%%%% Last updated: Nov 26 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear all; close all;
StartINTLAB();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PARAMETERS
%compare_SS_file = 'SS.mat';
%load(compare_SS_file);	%TMP PLACEMENT
%state_id = 1;
%ID_COLORS = {'b', 'b', 'm', 'm', 'r', 'r', 'y', 'y', 'k'};	%For the standard parameters (0.07, 0.025)

%THIS IS FOR THE LOCALIZED ATOMS STATE.
compare_SS_file = 'SS_LOCALIZED_FINAL/SS.mat';
state_id = 12;
id_color = 'r';
load(compare_SS_file);
list_SS{state_id}.A = list_SS{state_id}.A(1:41, 1:41);
pfc_g = PFCGeneralData(pfc_g.two_q, 40, pfc_g.Nx, 1, 16, pfc_g.intval_psibar, pfc_g.intval_beta, pfc_g.intval_nu);

%THIS IS FOR THE GRAIN BOUNDARY.
% compare_SS_file = 'Long_SS/SS_LargeDomain/SS.mat';
% id_color = 'r';
% load(compare_SS_file);
% list_SS{state_id}.A = list_SS{state_id}.A(1:41, 1:41);
% pfc_g = PFCGeneralData(pfc_g.two_q, 40, pfc_g.Nx, 1, 16, pfc_g.intval_psibar, pfc_g.intval_beta, pfc_g.intval_nu);


%Branch parameters
initial_psibar_direction = -1;
step_count = 900;

%Newton tolerance and adaptive arclength with a built-in maximum step size
norm_tolerance = 10e-14;	%fprintf('WARNING USING LOWER TOLERANCE FOR LOCALIZED'); %THIS SHOULD BE 10e-14 IN MOST CASES
max_s = 0.001;				%This should be lower when doing snakes

%Set to a NaN to disable showing the phase
mod_show_phase = 6;

%Don't change this sigh
delta_s = 0.001;
keller_formula = @(delta, iters) min(max_s, delta * 2^((4-iters)/3));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load the SS file and parameters
psibar_limit = sqrt(pfc_g.beta/1.5);
%id_color = ID_COLORS{state_id};

%Compute the current position and fix the initial "tangent" direction to the desired direction
cur_A = list_SS{state_id}.A;
M = length(cur_A) - 1;
cur_X = [pfc_g.psibar; SwitchMat2Vec(cur_A)];
Xdot = 0*cur_X; Xdot(1) = initial_psibar_direction;

%Compute dF, dpsibar and the initial derivatives
dFdpsibar = zeros((1+M)^2, 1); dFdpsibar(1) = -1;
G_PFC = DF_2D(cur_A, pfc_g, false);
G_X = [dFdpsibar, G_PFC];

%Accumulation arrays and initial values
psibar_arr = pfc_g.psibar; energy_arr = list_SS{state_id}.E; norm_arr = NormL2(cur_A);
A_arr = {};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Prepare the phase plot if needed
if ~isnan(mod_show_phase)
	[XX, YY] = meshgrid(pfc_g.x, pfc_g.y); plot_psi = 0*XX;
	jvec = 2*pi/pfc_g.Lx * (0:M); kvec = 2*pi/pfc_g.Ly * (0:M); weight = WeightMatrix(M);
	
	phase_fig = figure(1); phase_img = imagesc(pfc_g.x, pfc_g.y, plot_psi); colormap(flipud(bone));
	title(sprintf('Step %d, m = %.3f', 0, pfc_g.psibar));
	set(gca, 'XTickLabel', []); set(gca, 'YTickLabel', []); set(gca,'YDir','normal'); axis equal;
end

%Iterate over the branch
for iter = 2:step_count
	%Compute the tangent vector (only one is found!) and ensure we move in the same direction as before
	kernel = null(G_X);
	Xdot = kernel * sign(dot(kernel, Xdot));
	
	%Compute the predictor
	Xhat = cur_X + delta_s * Xdot; cur_X = Xhat;
	
	%Compute the corrector using the pseudo arclength Newton method
	newton_iters = 0;
	flag_bad_newton = true;
	while newton_iters < 200
		%Compute the matrices
		cur_A = SwitchVec2Mat(cur_X(2:end));
		pfc_g.psibar = cur_X(1);
		F_PAL = [(cur_X - Xhat)' * Xdot; F_2D(cur_A, pfc_g, false)];
		G_PFC = DF_2D(cur_A, pfc_g, false);
		G_X = [dFdpsibar, G_PFC];
		
		%Test the norm of the PAL F
		if norm(F_PAL, Inf) < norm_tolerance
			flag_bad_newton = false;
			break;
		end
		
		%Update the current X	
		DF_PAL = [Xdot'; G_X];		
		cur_X = cur_X - DF_PAL \ F_PAL;
		newton_iters = newton_iters + 1;
	end
	
	%Print a sad message
	if flag_bad_newton
		fprintf('OOPSIE! NEWTON DID NOT CONVERGE!!!');
	end
	
	%Adapt the arclength step
	delta_s = keller_formula(delta_s, newton_iters);
	
	%Accept the step along the branch
	cur_A = SwitchVec2Mat(cur_X(2:end));
	pfc_g.psibar = cur_X(1);
	G_PFC = DF_2D(cur_A, pfc_g, false);
	G_X = [dFdpsibar, G_PFC];
	
	%Compute the new energy, norm, determinant and Morse index
	psibar_arr(end+1) = pfc_g.psibar; energy_arr(end+1) = GetEnergy(cur_A,pfc_g,false); norm_arr(end+1) = NormL2(cur_A);
	
	%Save the array every step, debug only...
	if ~mod(iter, 5)
		A_arr{end+1} = cur_A;
	else
		A_arr{end+1} = [];
	end
	
	%Stop when psibar is too large
	if abs(psibar_arr(iter)) > psibar_limit break; end
	
	%Show the phase along the branch if needed, copied and slightly modified from GetPhase
	if mod(iter, mod_show_phase) == 0
		plot_psi = 0*plot_psi;
		for cur_j = 0:M
			for cur_k = 0:M
				plot_psi = plot_psi + weight(1+cur_k, 1+cur_j) * cur_A(1+cur_k, 1+cur_j) .* ...
					cos(jvec(1+cur_j) * XX) .* cos(kvec(1+cur_k) * YY);		
			end
		end
		set(phase_img, 'CData', plot_psi);
		phase_fig.Children.Title.String = sprintf('Step %d, Delta = %.4f, m = %.3f', iter, delta_s, pfc_g.psibar);
		pause(0.01);
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot the relative energy
figure(2); box on; hold on; xlabel('$\bar{\psi}$', 'Interpreter', 'latex'); ylabel('Energy');
plot([-psibar_limit, psibar_limit], [0, 0], 'k', 'linewidth', 2);
const_energy_arr = psibar_arr.^2/2 + (psibar_arr.^2 - pfc_g.beta).^2/4;
plot(psibar_arr, energy_arr - const_energy_arr, 'color', id_color, 'linewidth', 2.0);
%scatter(psibar_arr(1), energy_arr(1)-const_energy_arr(1), 20, 'k', 'fill');
%xlim([-0.15, 0.15]); ylim([-12e-5, 1e-5]);		%For standard parameters
xlim([0.4, 0.7]); ylim([-1e-2, 1e-2]);			%LOCALIZED

%Plot the norm
figure(3); box on; hold on; xlabel('$\bar{\psi}$', 'Interpreter', 'latex'); ylabel('Norm');
plot([-psibar_limit, 0, psibar_limit], [psibar_limit, 0, psibar_limit], 'k', 'linewidth', 2);
plot(psibar_arr, norm_arr, 'color', id_color, 'linewidth', 2.0);
%scatter(psibar_arr(1), norm_arr(1), 20, 'k', 'fill');
%xlim([-0.15, 0.1501]); ylim([0.0, 0.15]);		%For standard parameters
xlim([0.4, 0.65]); ylim([0.4, 0.6]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%To save images for publications, save plots with h(n) = plot(...)
%This must be done for each figure.

% legend(h, 'Constant', 'Atoms', 'Stripes', 'Donuts', 'Rectangles')
 set(findall(gcf, '-property', 'FontSize'), 'FontSize', 16)
 ax = gca; outerpos = ax.OuterPosition; ti = ax.TightInset;
 ax.Position = [outerpos(1)+ti(1), outerpos(2)+ti(2), outerpos(3)-ti(1)-ti(3)-0.01, outerpos(4)-ti(2)-ti(4)];
 set(gcf, 'Color', 'w'); set(gcf, 'InvertHardcopy', 'off');
% saveas(1, 'Continue_Energy.png');
