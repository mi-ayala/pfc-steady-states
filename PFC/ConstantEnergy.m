function E0 = ConstantEnergy(two_q, psibar, beta)
	%Constant energy
	if isnan(two_q)
		E0 = (1-beta)*psibar^2/2 + psibar^4/4 + beta^2/4;
	else
		E0 = (two_q^4-beta)*psibar^2/2 + psibar^4/4 + beta^2/4;
	end
end
