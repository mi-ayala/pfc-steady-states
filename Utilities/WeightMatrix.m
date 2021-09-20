%Weights accounting for the symmetries given by the Neumann boundary conditions
function weight = WeightMatrix(M)
	weight = 4 * ones(M+1);
	weight(:, 1) = 2;
	weight(1, :) = 2;
	weight(1, 1) = 1;
end
