%Sorts the steady states according to their energy
function Sort_SS(compare_SS_file)
	load(compare_SS_file);
	
	%Get the sort order
	e_extract = cellfun(@(obj)obj.E, list_SS);
	[~, nn] = sort(e_extract);
	if isequal(nn, 1:length(list_SS))
		return;
	end
	
	%Sort and save
	list_SS = list_SS(nn);
	list_counts = list_counts(nn);
	save(compare_SS_file, '-append', 'list_SS', 'list_counts');
	
	%Replot the phases
	for nn = 1:length(list_SS)
		figure('visible','off');
		SavePhaseImage(list_SS{nn}.psi, pfc_g, [fileparts(compare_SS_file), sprintf('/SS_%d.png', nn)]);
	end
end
