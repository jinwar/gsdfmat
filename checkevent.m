% This script can only be run after average_event.m is run.
%
% written by ge jin, 
% jinwar@gmail.com
% Mar, 2012

eventn=input('Event Num:');

	figure(5)
	clf
	hold on
		ax = usamap(lalim, lolim);
		set(ax, 'Visible', 'off')
		states = shaperead('usastatehi', 'UseGeoCoords', true);
		geoshow(ax, states, 'FaceColor', [0.7 0.7 0.7])
		surfacem(xi,yi,eventdata(eventn).GV')
%		plotm(stadata(:,2),stadata(:,3),'rv');
		%contourm(xi,yi,z,30,'k');
		geoshow(ax, states, 'FaceColor', 'none')
		seiscolormap
		load seiscmap
		colormap(seiscmap);
		caxis([avgphv*0.95 avgphv*1.05]);
		colorbar
		filename = sprintf('isomap_dyn_%s_%3d',eventdata(eventn).eventid,periods(period+1));
		title(filename,'Interpreter','none');

	figure(6)
	clf
	hold on
		ax = usamap(lalim, lolim);
		set(ax, 'Visible', 'off')
		states = shaperead('usastatehi', 'UseGeoCoords', true);
		geoshow(ax, states, 'FaceColor', [0.7 0.7 0.7])
		surfacem(xi,yi,eventdata(eventn).GV_correct')
%		plotm(stadata(:,2),stadata(:,3),'rv');
		%contourm(xi,yi,z,30,'k');
		geoshow(ax, states, 'FaceColor', 'none')
		seiscolormap
		load seiscmap
		colormap(seiscmap);
		caxis([avgphv*0.95 avgphv*1.05]);
		colorbar
		filename = sprintf('isomap_oricor_%s_%3d',eventdata(eventn).eventid,periods(period+1));
		title(filename,'Interpreter','none');

	figure(7)
	clf
	hold on
		ax = usamap(lalim, lolim);
		set(ax, 'Visible', 'off')
		states = shaperead('usastatehi', 'UseGeoCoords', true);
		geoshow(ax, states, 'FaceColor', [0.7 0.7 0.7])
		surfacem(xi,yi,eventdata(eventn).GV_bestcor')
%		plotm(stadata(:,2),stadata(:,3),'rv');
		%contourm(xi,yi,z,30,'k');
		geoshow(ax, states, 'FaceColor', 'none')
		seiscolormap
		load seiscmap
		colormap(seiscmap);
		caxis([avgphv*0.95 avgphv*1.05]);
		colorbar
		filename = sprintf('isomap_bestcor_%s_%3d',eventdata(eventn).eventid,periods(period+1));
		title(filename,'Interpreter','none');

	figure(8)
	clf;
		plot(alpha,eventdata(eventn).err);

	figure(9)
	clf
	hold on
		ax = usamap(lalim, lolim);
		set(ax, 'Visible', 'off')
		states = shaperead('usastatehi', 'UseGeoCoords', true);
		geoshow(ax, states, 'FaceColor', [0.7 0.7 0.7])
		surfacem(xi,yi,eventdata(eventn).GV_correct'-eventdata(eventn).GV')
%		plotm(stadata(:,2),stadata(:,3),'rv');
		%contourm(xi,yi,z,30,'k');
		geoshow(ax, states, 'FaceColor', 'none')
		seiscolormap
		load seiscmap
		colormap(seiscmap);
%		caxis([avgphv*0.95 avgphv*1.05]);
		colorbar
		filename = sprintf('oricor_%s_%3d',eventdata(eventn).eventid,periods(period+1));
		title(filename,'Interpreter','none');

	figure(10)
	clf
	hold on
		ax = usamap(lalim, lolim);
		set(ax, 'Visible', 'off')
		states = shaperead('usastatehi', 'UseGeoCoords', true);
		geoshow(ax, states, 'FaceColor', [0.7 0.7 0.7])
		surfacem(xi,yi,eventdata(eventn).GV_correct'-isoGV')
%		plotm(stadata(:,2),stadata(:,3),'rv');
		%contourm(xi,yi,z,30,'k');
		geoshow(ax, states, 'FaceColor', 'none')
		seiscolormap
		load seiscmap
		colormap(seiscmap);
%		caxis([avgphv*0.95 avgphv*1.05]);
		colorbar
		filename = sprintf('oricor_to_isomap_%s_%3d',eventdata(eventn).eventid,periods(period+1));
		title(filename,'Interpreter','none');

