%Computes the (averaged) dot product in L2 of two coefficient matrices
function dl2 = DotL2(A, B)
	dl2 = sum(sum(WeightMatrix(length(A)-1) .* A .* B));
end
