%Converts a matrix F into its reshaped vector formulation
function vec_F = SwitchMat2Vec(mat_F)
	M = length(mat_F)-1;
	vec_F = reshape(mat_F', [(1+M)^2, 1]);
end
