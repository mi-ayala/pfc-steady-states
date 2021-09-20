%Converts a vector F into its original matrix formulation
function mat_F = SwitchVec2Mat(vec_F)
	M = sqrt(length(vec_F))-1;
	mat_F = reshape(vec_F', [1+M, 1+M])';
end
