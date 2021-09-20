%Computes the signature of G and its unstable directions
function [num_pos_eig, num_neg_eig, num_zero_eig, unstable_directions, USE_INTLAB] = GetStability(G, USE_INTLAB)

if USE_INTLAB
	if ~isa(G, 'intval') error('WARNING!!! G is not an intval but USE_INTLAB is true!'); end
end

%Compute the eigenstuff
num_pos_eig = 0; num_neg_eig = 0; num_zero_eig = 0; unstable_directions = {};
[mat_eigenvectors, vec_eigenvalues] = eig(mid(G), 'vector');

%Only test eigenvalues that aren't too negative, this is necessary for large M, especially when using INTLAB
large_negative = vec_eigenvalues < -10;
if USE_INTLAB
	large_negative = vec_eigenvalues < -0.5;
end
mat_eigenvectors(:, large_negative) = []; vec_eigenvalues(large_negative) = [];

%Test each eigenvalue, using verifyeig when needed
validated_eigenvalues = false(size(vec_eigenvalues));
for nn = 1:length(vec_eigenvalues)
	fprintf('%d%% ', floor(100*nn/length(vec_eigenvalues)));
	
	if USE_INTLAB
		%Skip eigenvalues that have already been validated as part of a previous cluster
		if validated_eigenvalues(nn)
			continue;
		end
		
		%Identify clusters of equal eigenvalues
		cur_eigenvalue = vec_eigenvalues(nn);
		same_eigenvalues = [nn];
		inv_subspace = mat_eigenvectors(:, nn);
		for mm = (nn+1):length(vec_eigenvalues)
			if abs(cur_eigenvalue - vec_eigenvalues(mm)) < 10^-10
				same_eigenvalues(end+1) = mm;
				inv_subspace(:, end+1) = mat_eigenvectors(:, mm);
			end
		end
		
		%Verify the cluster of eigenvalues on the eigenspace
		[verified_eigenvalue, verified_eigenvectors] = verifyeig(G, cur_eigenvalue, inv_subspace);
		
		%Detect problems with verifyeig and give up. This sometimes occurs for diagonal matrices, in which case the
		%eigenvalues can be checked explicitly from the interval definition of G.
		if isnan(verified_eigenvalue)
			fprintf('Verifyeig output NaN, giving up.\n');
			num_pos_eig = NaN; num_neg_eig = NaN; num_zero_eig = NaN;
			USE_INTLAB = false;
			return;
		end
		
		%Check stability on the cluster
		for adopt_id = 1:length(same_eigenvalues)
			validated_eigenvalues(same_eigenvalues(adopt_id)) = true;
			
			%Check the sign of the eigenvalue(s)
			if real(verified_eigenvalue) > 0
				num_pos_eig = num_pos_eig + 1;
				unstable_directions{end+1} = SwitchVec2Mat(mid(verified_eigenvectors(:, adopt_id)));
			elseif real(verified_eigenvalue) < 0
				num_neg_eig = num_neg_eig + 1;			
			else
				num_zero_eig = num_zero_eig + 1;
			end
		end
	else
		%Simply check the sign of the eigenvalue
		if vec_eigenvalues(nn) > 0
			num_pos_eig = num_pos_eig + 1;
			unstable_directions{end+1} = SwitchVec2Mat([mat_eigenvectors(:, nn)]);
		elseif vec_eigenvalues(nn) < 0
			num_neg_eig = num_neg_eig + 1;			
		else
			num_zero_eig = num_zero_eig + 1;
		end
	end
end

%TEMPORARY VERIFICATION IN CASE I SCREWED UP WITH INTLAB - This still fails sometimes, why does verifyeig output NaN?
VERIF_num_pos_eig = 0; VERIF_num_neg_eig = 0; VERIF_num_zero_eig = 0;
for nn = 1:length(vec_eigenvalues)
	if vec_eigenvalues(nn) > 0
		VERIF_num_pos_eig = VERIF_num_pos_eig + 1;
	elseif vec_eigenvalues(nn) < 0
		VERIF_num_neg_eig = VERIF_num_neg_eig + 1;			
	else
		VERIF_num_zero_eig = VERIF_num_zero_eig + 1;
	end
end
assert(isequal([num_pos_eig, num_neg_eig, num_zero_eig], [VERIF_num_pos_eig, VERIF_num_neg_eig, VERIF_num_zero_eig]));


%Output verification state, there should be a single unstable direction for stable states
if USE_INTLAB fprintf('(VERIFIED) '); end

num_neg_eig = num_neg_eig + sum(large_negative);
fprintf('Skipped %d large negative eigenvalues.\n', sum(large_negative));

%Output the stability assessment
if num_zero_eig == 0
	if num_pos_eig == 1 fprintf('The steady state is stable (in the span of cosines).\n');
	else fprintf('There are %d unstable directions.\n', num_pos_eig-1); end
else
	fprintf('There are %d vanishing/small eigenvalues!\n', num_zero_eig);
end

end
