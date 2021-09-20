function PhaseDiagram_Show(xx, yy, show_pd, FILE_NAME)
	%Keep a copy of the original
	original_show_pd = show_pd;
	%xx = 1:28; yy = 1:28;
	%Replace NaN with the closest value, assuming we properly enclosed the transition regions and edges
	show_pd = [NaN(size(show_pd, 1), 1), show_pd, NaN(size(show_pd, 1), 1)];
	show_pd = [NaN(1, size(show_pd, 2)); show_pd; NaN(1, size(show_pd, 2))];
	[bad_x, bad_y] = find(isnan(show_pd(2:end-1, 2:end-1)));
	while ~isempty(bad_x)
		mod_show_pd = show_pd;
		for nn = 1:length(bad_x)
			neighbors = [show_pd(bad_x(nn)+2, bad_y(nn)+1), show_pd(bad_x(nn)+1, bad_y(nn)+2), ...
				show_pd(bad_x(nn), bad_y(nn)+1), show_pd(bad_x(nn)+1, bad_y(nn))];
			mod_show_pd(bad_x(nn)+1, bad_y(nn)+1) = mode(neighbors);
		end
		show_pd = mod_show_pd;
		[bad_x, bad_y] = find(isnan(show_pd(2:end-1, 2:end-1)));
	end
	show_pd = show_pd(2:end-1, 2:end-1);
	
	%Show all other modes as black
	original_show_pd(original_show_pd > 3) = 4;
	show_pd(show_pd > 3) = 4;
	
	%Plot phase diagram
	figure(1); hold on; box on; imagesc(xx, yy, show_pd, 'AlphaData', 0.5); set(gca, 'YDir', 'normal')
	xlim([min(xx), max(xx)]); ylim([min(yy), max(yy)]);
	xlabel('$\bar{\psi}$', 'Interpreter', 'latex'); ylabel('$\beta$', 'Interpreter', 'latex');
	cmap = [1,0,0; 0,0,1; 1,1,0; 0,1,0]; colormap(cmap(1:(max(show_pd(:))), :));
	
	%Points
	h = zeros(4, 1);
	[py, px] = find(original_show_pd == 1); h(1) = scatter(xx(px), yy(py), 25, 'r', 'filled', 'MarkerEdgeColor', 'k');
	[py, px] = find(original_show_pd == 2); h(2) = scatter(xx(px), yy(py), 25, 'b', 'filled', 'MarkerEdgeColor', 'k');
	[py, px] = find(original_show_pd == 3); h(3) = scatter(xx(px), yy(py), 25, 'y', 'filled', 'MarkerEdgeColor', 'k');
	[py, px] = find(original_show_pd == 4); h(4) = scatter(xx(px), yy(py), 25, 'g', 'filled', 'MarkerEdgeColor', 'k');
	
	%Curves
	finer_xx = linspace(min(xx), max(xx), 200);
	plot(finer_xx, 20.22*finer_xx.^2, 'k', 'linewidth', 1.5);
	plot(finer_xx, 37/15*finer_xx.^2, 'w', 'linewidth', 1.5);
	
	%Custom legend
	%h(5) = plot(NaN, NaN, 'w', 'linewidth', 2);
	%lgnd = legend(h, 'Constant', 'Stripes', 'Atoms', 'Others', '$\beta = 37/15 \bar{\psi}^2$');
	h(5) = plot(NaN, NaN, 'k', 'linewidth', 2); h(6) = plot(NaN, NaN, 'w', 'linewidth', 2);
	lgnd = legend(h, 'Constant', 'Stripes', 'Atoms', 'Others', '$\beta = 20.22\bar{\psi}^2$', '$\beta = 37/15 \bar{\psi}^2$');
	set(lgnd, 'Location', 'Northwest'); set(lgnd, 'Interpreter', 'latex'); set(lgnd, 'Color', 0.7*[1,1,1]);
	
	%Printing
	set(findall(gcf, '-property', 'FontSize'), 'FontSize', 16)
	ax = gca; outerpos = ax.OuterPosition; ti = ax.TightInset;
	ax.Position = [outerpos(1)+ti(1), outerpos(2)+ti(2), outerpos(3)-ti(1)-ti(3), outerpos(4)-ti(2)-ti(4)];
	set(gcf, 'Color', 'w'); set(gcf, 'InvertHardcopy', 'off');
	saveas(1, FILE_NAME);
end
