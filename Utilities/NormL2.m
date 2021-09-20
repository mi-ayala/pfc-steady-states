%Computes the (average) L2 norm of a coefficient matrix
function nl2 = NormL2(A)
	A = A.^2;
	A(:, 2:end) = 2*A(:, 2:end);
	A(2:end, :) = 2*A(2:end, :);
	nl2 = sqrt(sum(A(:)));
end
