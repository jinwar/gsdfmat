clear

% Input event list file
eventslist='testevent';

% some constants
ampmintor=0.5;
isonefigure=1;


eventfp=fopen(eventslist,'r');

event=fgetl(eventfp)

while ischar(event)
for period=0:0

filename=sprintf('%s_%1d.mat',event,period);
load(filename);

stlamid=mean(lalim);
stlomid=mean(lolim);

azi_array=azimuth(stlamid,stlomid,evla,evlo);

for i=1:m
	for j=1:n
		if xi(i,j) < 42 || xi(i,j)>48
			%dz_correct(i,j)=NaN;
		end
	end
end

lamda=100:10:300;
theta=azi_array-20:2:azi_array+20;

%cwt=CWT_geomorlet2d(xi,yi,dz_correct,185,azi);
for i=1:length(lamda)
	for j=1:length(theta)
		fcwt=FFTCWT_geomorlet2d(xi,yi,dz_correct,lamda(i),theta(j));
		cwtmaxamp(i,j)=max(max(abs(fcwt)));
	end
end

[mi mj]=find(cwtmaxamp==max(max(cwtmaxamp)));
%best_cwt=CWT_geomorlet2dV2(xi,yi,dz_correct,lamda(mi),aspect);
best_cwt=FFTCWT_geomorlet2d(xi,yi,dz_correct,lamda(mi),theta(mj));
[dz_cwtcor,amp]=cwtgeomorlet2d_correct(dz_correct,best_cwt);


[lamdax thetay]=meshgrid(lamda,theta);
figure(1)
clf

% 	subplot('position',[0.05 0.1 0.5 0.8])
	hold on
	contour(lamdax,thetay,cwtmaxamp',30);
	plot(lamda(mi),theta(mj),'x','Markersize',20)
	xlabel('Interference Wave Length/km','fontsize',20);
	ylabel('Azimuth','fontsize',20);
	stemp=sprintf('Wavelength: %5.1f km, Azimuth: %5.1f, amp: %f', lamda(mi), theta(mj),amp);
	h=title(stemp);
	set(h,'fontsize',20)
	shading flat;
	colorbar
    filename=sprintf('%s_%1d_grd_cwtsearch',event,period)
    print('-depsc',filename)
	

% The following is to plot the CWT result
figure(2)
clf
% 	subplot('position',[0.6 0.05 0.3 0.4])
	ax = usamap(lalim, lolim);
	set(ax, 'Visible', 'off')
	states = shaperead('usastatehi', 'UseGeoCoords', true);
	geoshow(ax, states, 'FaceColor', [0.5 0.5 1])
	surfacem(xi,yi,abs(best_cwt));
	%contourm(xi,yi,z,30,'k');
	geoshow(ax, states, 'FaceColor', 'none')
    colorbar
	seiscolormap
	load seiscmap
	colormap(seiscmap);
    filename=sprintf('%s_%1d_grd_cwtamp',event,period)
    print('-depsc',filename)


% Plot the best fcwt
figure(3)
clf
%     subplot('position',[0.6 0.5 0.3 0.4])
	ax = usamap(lalim, lolim);
	set(ax, 'Visible', 'off')
	states = shaperead('usastatehi', 'UseGeoCoords', true);
	geoshow(ax, states, 'FaceColor', [0.5 0.5 1])
	surfacem(xi,yi,dz_cwtcor);
	colorbar
	seiscolormap
	load seiscmap
	colormap(seiscmap);
	caxis(phvrange(period+1,:));
    filename=sprintf('%s_%1d_grd_cwtcor',event,period)
    print('-depsc',filename)

figure(4)
clf
% 	subplot('position',[0.6 0.05 0.3 0.4])
	ax = usamap(lalim, lolim);
	set(ax, 'Visible', 'off')
	states = shaperead('usastatehi', 'UseGeoCoords', true);
	geoshow(ax, states, 'FaceColor', [0.5 0.5 1])
	surfacem(xi,yi,dz_correct);
	contourm(xi,yi,z,30,'k');
	geoshow(ax, states, 'FaceColor', 'none')
	colorbar
	seiscolormap
	load seiscmap
	colormap(seiscmap);
	caxis(phvrange(period+1,:));
    filename=sprintf('%s_%1d_grd_ampcor',event,period)
    print('-depsc',filename)

% figure(1)
%     filename=sprintf('cwtplot_%s_%1d',event,period);
% 	%print('-dpng',filename);
% 	export_fig(filename,'-jpg');


end   % end of periods

event=fgetl(eventfp)

end	  % end of events
