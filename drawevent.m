function y=drawevent(eventname, period)

drawcmt=0;
filename=sprintf('%s_%1d.mat',eventname,period);
load(filename);

figure(1);
clf
	subplot('position',[0.45 0.55 0.15 0.4]);
	ax = usamap(lalim, lolim);
	set(ax, 'Visible', 'off')
	states = shaperead('usastatehi', 'UseGeoCoords', true);
	geoshow(ax, states, 'FaceColor', [0.5 0.5 1])
	surfacem(xi,yi,dz);
	contourm(xi,yi,z,30,'k');
	geoshow(ax, states, 'FaceColor', 'none')
	seiscolormap
	load seiscmap
	colormap(seiscmap);
	%quiverm(bigxi,bigyi,u,v);
	caxis(phvrange(period+1,:));
	colorbar;
	%plotm(coor(:,1),coor(:,2),'rv');
	filename = sprintf('%s_%1d_grdmap',event,period);
	title(filename,'Interpreter','none');

% Calculate the aspect change compare to the event back azimuth
%figure(2);
%clf
	subplot('position',[0.05 0.05 0.15 0.4]);
    ax = usamap(lalim, lolim);
    set(ax, 'Visible', 'off')
    states = shaperead('usastatehi', 'UseGeoCoords', true);
    geoshow(ax, states, 'FaceColor', [0.5 0.5 1])
	surfacem(xi,yi,aspectch);
	contourm(xi,yi,z,30,'k');
	geoshow(ax, states, 'FaceColor', 'none')
	sc=nanmean(nanstd(aspectch));
	caxis([-4*sc 4*sc]);
	colorbar;
	plotm(coor(:,1),coor(:,2),'rv');
    filename=sprintf('Surface Surface Gradient Direction event:%s PID:%d',event,period);
    title(filename,'Interpreter','none');
	filename = sprintf('%s_%1d_azimap',event,period);
	%print('-depsc',filename);

%figure(3);
%clf
	subplot('position',[0.25 0.05 0.15 0.4]);
    ax = usamap(lalim, lolim);
    set(ax, 'Visible', 'off')
    states = shaperead('usastatehi', 'UseGeoCoords', true);
    geoshow(ax, states, 'FaceColor', [0.5 0.5 1])
	surfacem(xi,yi,ampmap);
	%caxis([2 5]);
	colorbar;
	%plotm(coor(:,1),coor(:,2),'rv');
    filename=sprintf('Amplitude Map event:%s PID:%d',event,period);
    title(filename,'Interpreter','none');
	filename = sprintf('%s_%1d_ampmap',event,period);
	%print('-depsc',filename);

%figure(4)
%clf
	subplot('position',[0.45 0.05 0.15 0.4]);
    ax = usamap(lalim, lolim);
    set(ax, 'Visible', 'off')
    states = shaperead('usastatehi', 'UseGeoCoords', true);
	surfacem(xi,yi,correctV);
	%contourm(xi,yi,z,30,'k');
	geoshow(ax, states, 'FaceColor', 'none')
	sc=nanmean(nanstd(correctV));
	caxis([-4*sc 4*sc]);
	seiscolormap
	load seiscmap
	colormap(seiscmap);
	colorbar;
	%plotm(coor(:,1),coor(:,2),'rv');
    filename=sprintf('Amplitude Correction Map event:%s PID:%d',event,period);
    title(filename,'Interpreter','none');
	filename = sprintf('%s_%1d_ampcorrection',event,period);
	%print('-depsc',filename);

%figure(5);
%clf
	subplot('position',[0.65 0.1 0.3 0.8])
	ax = usamap(lalim, lolim);
	set(ax, 'Visible', 'off')
	states = shaperead('usastatehi', 'UseGeoCoords', true);
	geoshow(ax, states, 'FaceColor', [0.5 0.5 1])
	surfacem(xi,yi,dz_correct);
	contourm(xi,yi,z,20,'k');
	geoshow(ax, states, 'FaceColor', 'none')
	seiscolormap
	load seiscmap
	colormap(seiscmap);
	%quiverm(bigxi,bigyi,u,v);
	caxis(phvrange(period+1,:));
	colorbar;
	plotm(coor(:,1),coor(:,2),'bv','markersize',5);
	filename = sprintf('%s_%1d_grdmap_correct',event,period);
	title(filename,'Interpreter','none');
	%print('-depsc',filename);
    
%figure(6);
%clf

	subplot('position',[0.05 0.55 0.35 0.4])
    %define lat/long limits
	%
    latmin=min([evla coor(:,1)'])-10; % Bottom Latitude
    latmax=max([evla coor(:,1)'])+10;  % Top Latitude

    longmin=min([evlo coor(:,2)'-360])-10; % Left Longitude
    longmax=max([evlo coor(:,2)'-360])+10; % Right Logitude
    
    if (longmax-longmin)>180
        temp=longmax;
        longmax=longmin+360+20;
        longmin=temp-20;
    end

    Pline=10; % Grid spacing for latidue lines
    Mline=10; % Grid spacing for longitude lines

    %define event location
    elat=evla; % Event Latitude
    elong=evlo; % Event Longitude

    hh=axesm('mapproj','aitoff',...
       'maplatlim',[latmin latmax],'maplonlim',[longmin longmax],...
       'MLineLocation',Mline,'PLineLocation',Pline,'Grid','on',...
       'MeridianLabel','on','ParallelLabel','on');
    axis off, gridm on, framem on;
    geoshow('landareas.shp')
    
	if drawcmt
		cmtsolution = geteventcmt(eventname,evla, evlo)
	end

    hold on

    for i=1:length(coor)
       [la, lo]=gcwaypts(elat,elong,coor(i,1),coor(i,2),30);
     %  geoshow(slat(i),slong(i),'color','k','marker','.');
	   [dist azi]=distance(elat, elong, coor(i,1), coor(i,2));
	   if drawcmt
		   if isonnode(azi,cmtsolution)
			   hh=geoshow(la,lo,'displaytype','line','color','r');
		   else
			   hh=geoshow(la,lo,'displaytype','line','color','b');
		   end
	   else
		   hh=geoshow(la,lo,'displaytype','line','color','b');
	   end
       set(hh,'LineWidth',1)
    end
    
    [dist baz]=distance(40,-105,evla,evlo);
    filename=sprintf('Dist=%f, Baz=%f',dist,baz);
    title(filename);

% plot the fft2 spectrum
%	subplot('position',[0.15 0.55 0.25 0.4])
%	abfftdz=abs(fft2(nan2zero(dz_correct)));
%	surface(abfftdz);
%	shading flat;
%	caxis([0 mean(mean(abfftdz))]);
%	[m n]=size(abfftdz);
%	xlim([0 n])
%	ylim([0 m])
%	colorbar
	

figure(1)
    filename=sprintf('eventplot_%s_%1d',eventname,period);
	print('-dpng',filename);
	

