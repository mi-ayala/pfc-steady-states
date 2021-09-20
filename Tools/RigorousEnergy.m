%Computes the energy and error of the verified steady state near A
function [energy, energy_error] = RigorousEnergy(A, radius, pfc_g)
	%Compute the energy of the steady state with INTLAB
	A = intval(A);
	energy = GetEnergy(A, pfc_g, true);
	
	%Maybe implement this later
	if ~isnan(pfc_g.two_q)
		printf('Multimode PFC rigorous energy is not yet supported.\n');
		energy = mid(energy);
		energy_error = NaN;
		return
	end
	
	%Norm of A and its derivatives
	normA = NormMatrixNu(A, pfc_g.intval_nu_mat);
	normLaplA = NormMatrixNu(pfc_g.intval_lapl .* A, pfc_g.intval_nu_mat);
	normLaplLaplA = NormMatrixNu(pfc_g.intval_lapl.^2 .* A, pfc_g.intval_nu_mat);
	
	%Compute the sums
	rho = 1/pfc_g.intval_nu^2;
	S1 = abs(pfc_g.intval_lapl(2, 2)) * (rho^2+rho)/(1-rho)^4;
	S2 = ((pfc_g.intval_lapl(2, 1)^2 + pfc_g.intval_lapl(1, 2)^2) * (rho^4 + 11*rho^3 + 11*rho^2 + rho) + ...
		2*abs(pfc_g.intval_lapl(2, 1)*pfc_g.intval_lapl(1,2)) * (rho^4 + 2*rho^3 + rho^2))/(1-rho)^6;
	
	%Compute the energy radius
	energy_distance = (abs(1-pfc_g.intval_beta)*normA + normA^3 + 2*normLaplA + normLaplLaplA) * radius + ...
		0.5 * (2*S1 + S2 + abs(1-pfc_g.intval_beta) + 3*normA^2) *radius^2 + normA * radius^3 + 0.25 * radius^4;
	
	%Decouple the interval into the value and its radius
	energy_error = rad(energy) + sup(energy_distance);
	energy = mid(energy);
end
