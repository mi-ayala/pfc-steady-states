%Computes an approximate zero of F starting at the input A
function [A, failed_to_converge, iteration_step] = NewtonSolver_2D(A, pfc_g, max_iterations, tolerance, VERBOSE)
	if VERBOSE fprintf('Newtown solver, current norm: '); end
	
	%Start the iteration
	iteration_step = 0;
	failed_to_converge = true;
	vec_F = F_2D(A, pfc_g, false);
	cur_norm = sum(sum(abs(vec_F)));
	while iteration_step < max_iterations
		%Stop iterating if the sum of the coefficients of F are smaller than the threshold
		if VERBOSE fprintf('%.1e, ', cur_norm); end
		if cur_norm < tolerance
			failed_to_converge = false;
			break;
		end
		
		%Iterate according to the rule X -> X - (grad F)^-1 F, ignoring the first coefficient
		vec_F = F_2D(A, pfc_g, false);
		%temp_DF = DF_2D(A, pfc_g, false);
		%fprintf('Det: %.2e, (min, max) : (%.2e, %.2e), normF: %.2e\n', det(temp_DF), min(temp_DF(:)), max(temp_DF(:)), cur_norm);
		A = A - SwitchVec2Mat(DF_2D(A, pfc_g, false) \ vec_F);
		cur_norm = sum(sum(abs(vec_F)));
		iteration_step = iteration_step + 1;
	end
	
	if VERBOSE fprintf('\n'); end
end

