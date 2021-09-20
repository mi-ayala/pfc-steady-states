%Starts INTLAB if it isn't already
function StartINTLAB()
	try intval(0);
	catch
		fprintf('Starting Intlab.\n');
		cur_dir = cd; cd('~/Intlab/Intlab_V9/'); startintlab; cd(cur_dir); intvalinit('DisplayMidrad');
	end
end
