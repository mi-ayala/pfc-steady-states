%Obtain the Cosine Fourier coefficients from a phase field psi to a given truncation order M
%Also optionally computes the norm of the L2 components that cannot be represented by the cosines
function [A, ortho_psi] = GetCoefficients(psi, as_M, pfc_g, compute_orthonormal)
	if as_M ~= pfc_g.M printf('Output truncation order is different from original truncation.\n'); end
	
	%Fourier wavelengths
	[XX, YY] = meshgrid(pfc_g.x, pfc_g.y);
	jvec = 2*pi/pfc_g.Lx * (0:as_M);
	kvec = 2*pi/pfc_g.Ly * (0:as_M);
	
	%Integrate all the modes
	A = zeros(as_M+1);
	if ~compute_orthonormal
		for cur_j = 0:as_M
			for cur_k = 0:as_M	
				A(1+cur_k, 1+cur_j) = mean(mean(psi .* cos(jvec(1+cur_j) * XX) .* cos(kvec(1+cur_k) * YY)));
			end
		end
	else
		%Optionally also computes the norm of the orthogonal components if needed
		weight = WeightMatrix(as_M);
		ortho_psi = psi;
		for cur_j = 0:as_M
			for cur_k = 0:as_M	
				A(1+cur_k, 1+cur_j) = mean(mean(psi .* cos(jvec(1+cur_j) * XX) .* cos(kvec(1+cur_k) * YY)));
				ortho_psi=ortho_psi-weight(1+cur_k,1+cur_j)*A(1+cur_k,1+cur_j)*cos(jvec(1+cur_j)*XX).*cos(kvec(1+cur_k)*YY);	
			end
		end
	end
end
