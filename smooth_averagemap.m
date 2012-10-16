%% This script is used to smooth the output of average_event.m, and plot the smoothed maps
% Written by Ge Jin, jinwar@gmail.com
clear;

r0=[[-0.08 0.06]];

for iperiod=4:4
	
	disp('loading data');
	filename = sprintf('averageevent_%1d.mat',iperiod);
	load(filename);
    r=r0*(0.96^iperiod)
	disp('smoothing maps');
	D = max([100 2*periods(iperiod+1)]);
	sm_isoGV = smoothmap(xi,yi,newisoGV',D );
	sm_isoGV_bestcor = smoothmap(xi,yi,isoGV_bestcor', D);
	sm_cordiff = sm_isoGV_bestcor - sm_isoGV;
	nanid = mapnanid(xi,yi,newisoGV',D);
	sm_isoGV(nanid) = NaN;
	sm_isoGV_bestcor(nanid) = NaN;
	sm_cordiff(nanid) = NaN;

% 	disp('calculating std');
% 	for i=1:Nx
% 		disp(i/Nx)
% 		for j = 1:Ny
% 			clear phv phv_cor;
% 			for ie = 1:length(eventdata)
% 				phv(ie)=eventdata(ie).GV(i,j);
% 				phv_cor(ie) = eventdata(ie).GV_bestcor(i,j);
% 				if isnan(phv_cor(ie))
% 					phv(ie)=NaN;
% 				end
% 			end
% 			GVstd(i,j) = nanstd(phv);
% 			GVstd_bestcor(i,j) = nanstd(phv_cor);
% 		end
% 	end

	disp('drawing maps');
    figure(1)
    clf
    hold on
    ax = usamap(lalim, lolim);
    set(ax, 'Visible', 'off')
    states = shaperead('usastatehi', 'UseGeoCoords', true);
    geoshow(ax, states, 'FaceColor', [0.7 0.7 0.7])
    surfacem(xi,yi,sm_isoGV)
    geoshow(ax, states, 'FaceColor', 'none')
    load seiscmap
    colormap(seiscmap);
    colorbar
    avgphv=nanmean(nanmean(newisoGV));
    caxis([avgphv*(1+r)]);
    filename = sprintf('isomap_grd_%3d',periods(period+1));
    title(filename,'Interpreter','none');
    filename = sprintf('sm_avg_isomap_grd_%1d',period);
    print('-dpng',filename);
    
    figure(2)
    clf
    hold on
    ax = usamap(lalim, lolim);
    set(ax, 'Visible', 'off')
    states = shaperead('usastatehi', 'UseGeoCoords', true);
    geoshow(ax, states, 'FaceColor', [0.7 0.7 0.7])
    surfacem(xi,yi,sm_isoGV_bestcor);
    %		plotm(stadata(:,2),stadata(:,3),'rv');
    %contourm(xi,yi,z,30,'k');
    geoshow(ax, states, 'FaceColor', 'none')
    %seiscolormap
    load seiscmap
    colormap(seiscmap);
    caxis([avgphv*(1+r)]);
    colorbar
    filename = sprintf('isomap_cor_%3d',periods(period+1));
    title(filename,'Interpreter','none');
    filename = sprintf('sm_avg_isomap_bestcor_%1d',period);
    print('-dpng',filename);

    figure(3)
    clf
    hold on
    ax = usamap(lalim, lolim);
    set(ax, 'Visible', 'off')
    states = shaperead('usastatehi', 'UseGeoCoords', true);
    geoshow(ax, states, 'FaceColor', [0.7 0.7 0.7])
    surfacem(xi,yi,sm_cordiff)
    %		plotm(stadata(:,2),stadata(:,3),'rv');
    %contourm(xi,yi,z,30,'k');
    geoshow(ax, states, 'FaceColor', 'none')
    %seiscolormap
    load seiscmap
    colormap(seiscmap);
	avgdiff = nanmean(abs(sm_cordiff(:)));
    caxis([-avgdiff*2 avgdiff*2])
    colorbar
    filename = sprintf('cordiff_%3d',periods(period+1));
    title(filename,'Interpreter','none');
    filename = sprintf('sm_cordiff_%1d',period);
    print('-dpng',filename);

%     figure(4)
%     clf
%     hold on
%     ax = usamap(lalim, lolim);
%     set(ax, 'Visible', 'off')
%     states = shaperead('usastatehi', 'UseGeoCoords', true);
%     geoshow(ax, states, 'FaceColor', [0.7 0.7 0.7])
%     surfacem(xi,yi,GVstd'-GVstd_bestcor')
%     geoshow(ax, states, 'FaceColor', 'none')
% 	avgdiff = nanmean(abs(GVstd(:)));
% %    caxis([0 avgdiff*1])
%     colorbar
%     filename = sprintf('GVstd_%3d',periods(period+1));
%     title(filename,'Interpreter','none');
%     filename = sprintf('StdDiff_%1d',period);
%     print('-dpng',filename);
% 
%     figure(5)
%     clf
%     hold on
%     ax = usamap(lalim, lolim);
%     set(ax, 'Visible', 'off')
%     states = shaperead('usastatehi', 'UseGeoCoords', true);
%     geoshow(ax, states, 'FaceColor', [0.7 0.7 0.7])
%     surfacem(xi,yi,GVstd_bestcor')
%     %		plotm(stadata(:,2),stadata(:,3),'rv');
%     %contourm(xi,yi,z,30,'k');
%     geoshow(ax, states, 'FaceColor', 'none')
%     %seiscolormap
% 	avgdiff = nanmean(abs(GVstd(:)));
%     caxis([0 avgdiff*1])
%     colorbar
%     filename = sprintf('GVstd_bestcor_%3d',periods(period+1));
%     title(filename,'Interpreter','none');
%     filename = sprintf('GVstd_bestcor_%1d',period);
%     print('-dpng',filename);
end

	
