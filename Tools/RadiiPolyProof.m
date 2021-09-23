%Radii polynomial proof for the PFC equation
function [r_min, r_max, G] = RadiiPolyProof(A, pfc_g)
r_min = NaN; r_max = NaN;

%Compute the padded intval Fourier matrices
M = pfc_g.M;
[~, intval_gamma_pad, intval_nu_pad] = FourierMatrices(3*M, pfc_g, true);

%Solution norm and F (matrices)
A = intval(A);
norm_A = NormMatrixNu(A, pfc_g.intval_nu_mat);
FA = F_2D(A, pfc_g, true);

%DF and its inverse (tensors) - note that G^-1 will use interval arithmetic, but this is not necessary!
G = DF_2D(A, pfc_g, true);
Ginv = intval(inv(mid(G)));
if any(isnan(Ginv))
	fprintf('WARNING!!! G is not invertible!\n');
	return;
end

%Approximate inverse times L, with the maximum value of gamma after truncation
max_Gamma = max(1/(intval_gamma_pad(1, M+2)), 1/(intval_gamma_pad(M+2, 1)));
GinvLaplacian = Ginv * diag(SwitchMat2Vec(pfc_g.intval_lapl));


%--------------BOUNDS-----------------------------------------------------------
%Bound Z0
Z0 = NormTensorNu(eye((1+M)^2) - Ginv * G, pfc_g.intval_nu_mat);

%Bound Y0
GinvFA = ConvAk_2D(A, 3) ./ intval_gamma_pad;
GinvFA(1:1+M, 1:1+M) = SwitchVec2Mat(Ginv * FA);
Y0 = NormMatrixNu(GinvFA, intval_nu_pad);

%Bound Z2(r) = Z2_r0 + Z2_r1*r
norm_GinvLaplacian = NormTensorNu(GinvLaplacian, pfc_g.intval_nu_mat) + max_Gamma;
Z2_r0 = 6*norm_GinvLaplacian * norm_A;
Z2_r1 = 3*norm_GinvLaplacian;

%Bound Z1
AA = abs(ConvAk_2D(A, 2));
AA = [AA, zeros(1+2*M, 2*M)];
AA = [AA; zeros(2*M, 1+4*M)];

%Indices on {0...3M}^2 \ {0...M}^2 then on the full {-3M...3M}^2 \ {-M...M}^2. Use scatter to visualize [quad_1, quad_2]
[mesh_ff1, mesh_ff2] = meshgrid((M+1):(3*M), 0:M);		%This is the "top-right corner" of {0...3M}^2 \ {0...M}^2
[mesh_gg1, mesh_gg2] = meshgrid(0:(3*M), (M+1):(3*M));	%This is the "bottom-left" and "bottom-right" corners

quad_1 = [mesh_ff1(:); mesh_gg1(:)]; quad_2 = [mesh_ff2(:); mesh_gg2(:)];	%These form {0...3M}^2 \ {0...M}^2
quad_1 = [quad_1; quad_1; -quad_1; -quad_1]; quad_2 = [quad_2; -quad_2; quad_2; -quad_2];	%And this tiles it

quad_nu = [pfc_g.intval_nu.^(mesh_ff1(:)+mesh_ff2(:)); pfc_g.intval_nu.^(mesh_gg1(:)+mesh_gg2(:))];
quad_nu = repmat(quad_nu, [4,1]);

%Computing the phis
Z1_phi = 0*A;
for cur_tau = 0:M
	for cur_sigma = 0:M
		Z1_phi(1+cur_tau, 1+cur_sigma) = max(AA(sub2ind(size(AA), ...
			1+abs(cur_tau - quad_1), 1+abs(cur_sigma - quad_2))) ./ quad_nu);	   
	end
end

%And finally putting everything together
Z1 = 3*(NormMatrixNu(SwitchVec2Mat(GinvLaplacian * SwitchMat2Vec(Z1_phi)), pfc_g.intval_nu_mat) + ...
	max_Gamma * norm_A^2);


%--------------RADII POLYNOMIAL CHECK-------------------------------------------
%Finding the roots of the radii polynomial with INTLAB
intval_rad_poly = polynom([Z2_r1, Z2_r0, Z0+Z1-1, Y0]);
rad_poly = mid([Z2_r1, Z2_r0, Z0+Z1-1, Y0]);

%Iterate over the positive real roots of the polynomial and verify them
r_pos = intval([]);
for poly_roots = roots(rad_poly)'
	if isreal(poly_roots) && (poly_roots == 0 || poly_roots > 0)	%Sometimes the root is exactly 0 for the constant state
		r_pos(end+1) = verifypoly(intval_rad_poly, poly_roots);
	end
end

%If the approach works, we find exactly two positive roots with disjoint interval (if not, we don't know where they are!)
if length(r_pos) == 2 && isnan(intersect(r_pos(1), r_pos(2)))
	%Take the most stringent roots, i.e. the ones that produce the smallest interval of negativity
	r_min = min(sup(r_pos)); r_max = max(inf(r_pos));
	fprintf('Radii polynomial (%.1e, %.1e, %.1e, %.1e): r_* = %.2e and r^* = %.2e.\n', rad_poly, r_min, r_max);
else
	fprintf('WARNING!!! Radii polynomial (%.1e, %.1e, %.1e, %.1e), nu = %.2f failed!\n', rad_poly, pfc_g.nu);
end

end
