%Computes the k-ic discrete convolution A*A*...*A, autodetects if INTLAB should be used
function convAk = ConvAk_2D(A, k)
	M = length(A)-1;
	
	%Pad each dimension of A to the next power of 2 of the minimum padding (k-1)M+1 in the (+,+) quadrant
	len_pad = 2^(nextpow2(2*k*M+1)-1) - M;
	A_pad = [zeros(1+M, len_pad), fliplr(A(:, 2:end)), A, zeros(1+M, len_pad-1)];
	A_pad = [zeros(len_pad, size(A_pad, 2)); flipud(A_pad(2:end, :)); A_pad; zeros(len_pad-1, size(A_pad, 2))];
	A_pad = ifftshift(A_pad);
		
	if ~isa(A, 'intval')
		%Without interval arithmetic, we can directly use Matlab's fft and ifft functions
		convAk = ifft2(fft2(A_pad).^k);
	else
		%With interval arithmetic, verifyfft must be used instead
		convAk = verifyfft(transpose(verifyfft(verifyfft(transpose(verifyfft(A_pad, 1)), 1).^k, -1)), -1);
	end
	
	%Only return the first (1+kM)^2 non-zero components
	convAk = real(convAk(1:(1+k*M), 1:(1+k*M)));
end
