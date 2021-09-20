%Computes the discrete convolution A*B, autodetects if INTLAB should be used
function convAB = ConvAB_2D(A, B)
	M = length(A)-1;
	
	%Pad each dimension of A and B to the next power of 2 of the minimum padding (k-1)M+1 in the (+,+) quadrant
	len_pad = 2^(nextpow2(4*M+1)-1) - M;
	A_pad = [zeros(1+M, len_pad), fliplr(A(:, 2:end)), A, zeros(1+M, len_pad-1)];
	A_pad = [zeros(len_pad, size(A_pad, 2)); flipud(A_pad(2:end, :)); A_pad; zeros(len_pad-1, size(A_pad, 2))];
	A_pad = ifftshift(A_pad);
	B_pad = [zeros(1+M, len_pad), fliplr(B(:, 2:end)), B, zeros(1+M, len_pad-1)];
	B_pad = [zeros(len_pad, size(B_pad, 2)); flipud(B_pad(2:end, :)); B_pad; zeros(len_pad-1, size(B_pad, 2))];
	B_pad = ifftshift(B_pad);
	
	if ~isa(A, 'intval')
		%Without interval arithmetic, we can directly use Matlab's fft and ifft functions
		convAB = ifft2(fft2(A_pad) .* fft2(B_pad));
	else
		%With interval arithmetic, verifyfft must be used instead
		convAB = verifyfft(transpose(verifyfft(verifyfft(transpose(verifyfft(A_pad, 1)), 1) .* ...
			verifyfft(transpose(verifyfft(B_pad, 1)), 1), -1)), -1);
	end
	
	%Only return the first (1+M)^2 components
	convAB = real(convAB(1:(1+M), 1:(1+M)));
end
