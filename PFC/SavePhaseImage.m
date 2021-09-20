%Plots the phase psi using imagesc or imagewrite
function SavePhaseImage(psi, pfc_g, save_path)
	figure('visible','off');
	imagesc(pfc_g.x, pfc_g.y, psi); colormap(flipud(bone));
	set(gca,'xtick',[]); set(gca,'ytick',[]);
	set(gca,'YDir','normal'); axis equal;
	set(gca, 'visible', 'off');
	
	%Remove white space and save
	ax = gca; outerpos = ax.OuterPosition; ti = ax.TightInset;
	ax.Position = [outerpos(1)+ti(1), outerpos(2)+ti(2), outerpos(3)-ti(1)-ti(3), outerpos(4)-ti(2)-ti(4)];
	saveas(gca, save_path);
end
