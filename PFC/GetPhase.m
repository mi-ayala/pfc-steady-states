%Constructs psi and its Laplacian from A assuming Neumann boundary conditions
function psi = GetPhase(A, pfc_g, show_phase)
	as_M = length(A)-1;
	if as_M ~= pfc_g.M fprintf('Input truncation order is different from original truncation.\n'); end
	
	%Fourier wavelengths
	[XX, YY] = meshgrid(pfc_g.x, pfc_g.y);
	jvec = 2*pi/pfc_g.Lx * (0:as_M);
	kvec = 2*pi/pfc_g.Ly * (0:as_M);
	
	%Add up all of the modes
	psi = zeros(length(pfc_g.y), length(pfc_g.x));
	weight = WeightMatrix(as_M);
	for cur_j = 0:as_M
		for cur_k = 0:as_M	
			psi = psi + weight(1+cur_k, 1+cur_j) * A(1+cur_k, 1+cur_j) .* ...
				cos(jvec(1+cur_j) * XX) .* cos(kvec(1+cur_k) * YY);		
		end
	end
	
	if nargin > 2 && show_phase
		imagesc(pfc_g.x, pfc_g.y, psi); colormap(flipud(bone));
	end
end
