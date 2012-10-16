clear
eventslist='testevent';

fp=fopen(eventslist,'r');

event=fgetl(fp);
while ischar(event)
for period=1:1

clear coor data ptxyz xyz xi yi z interpz;

filename = sprintf('%s_%1d.mat',event,period);

load(filename);

ptxyz(:,1)=coor(:,2)-360;
ptxyz(:,2)=coor(:,1);
ptxyz(:,3)=coor(:,3);
save 'phasetemp.xyz' ptxyz -ASCII
system('csh phasemap.gmt');

lalim=[27 50.5];
lolim=[-114 -95];

lolim=lolim+360;

xyz=load('gmtphsurf.xyz');
xyz(:,1)=xyz(:,1)+360;
[xi,yi]=meshgrid(lalim(1):0.1:lalim(2), lolim(1):0.1:lolim(2)); 
z=griddata(xyz(:,2),xyz(:,1),xyz(:,3),xi,yi,'cubic');
interpz=griddata(coor(:,1),coor(:,2),coor(:,3),xi,yi,'cubic');

%figure(1);
%clf
%	usamap(lalim,lolim);
%	[topoZ, refvec] = etopo('ETOPO5.DAT',1,lalim,lolim);
%	geoshow(topoZ, refvec,'DisplayType','texturemap');
%	colormap(demcmap(topoZ));
%	plotm(coor(:,1),coor(:,2),'rv');
%	title('Topography');
%	filename = sprintf('%s_topo',event);
%	print('-dpng',filename);

figure(2);
clf
    ax = usamap(lalim, lolim);
    set(ax, 'Visible', 'off')
    states = shaperead('usastatehi', 'UseGeoCoords', true);
    geoshow(ax, states, 'FaceColor', [0.5 0.5 1])
	contourm(xi,yi,z,30);
	colorbar;
	plotm(coor(:,1),coor(:,2),'rv');
	title('Phase Surface');
	filename = sprintf('%s_%1d_phmap',event,period);
	print('-dpng',filename);

% Calculate the gradient map
[aspect, slope, gradN, gradE]=gradientm(xi,yi,z);
% Get rid of expolaration area
[m,n]=size(interpz);
for i=1:m
    for j=1:n
        if isnan(interpz(i,j))
%            z(i,j)=nan;
%            aspect(i,j)=nan;
%            slope(i,j)=nan;
%            gradN(i,j)=nan;
%            gradE(i,j)=nan;
        end
    end
end
% change the unit from slope to velocity
dz=tand(slope).^-1/1e3;
[bigxi,bigyi]=meshgrid(lalim(1):0.5:lalim(2), lolim(1):0.5:lolim(2)); 
u=interp2(xi,yi,gradN,bigxi,bigyi);
v=interp2(xi,yi,gradE,bigxi,bigyi);
%dz=tand(slope);

figure(3);
clf
    ax = usamap(lalim, lolim);
    set(ax, 'Visible', 'off')
    states = shaperead('usastatehi', 'UseGeoCoords', true);
    geoshow(ax, states, 'FaceColor', [0.5 0.5 1])
	surfacem(xi,yi,dz);
	%quiverm(bigxi,bigyi,u,v);
	%caxis([2 5]);
	colorbar;
	plotm(coor(:,1),coor(:,2),'rv');
	seiscolormap
	load seiscmap
	colormap(seiscmap);
	%quiverm(bigxi,bigyi,u,v);
	caxis(phvrange(period+1,:));
	colorbar;

	title('Surface Surface Gradient');
	filename = sprintf('%s_%1d_grdmap',event,period);
	%print('-dpng',filename);

figure(4);
clf
    ax = usamap(lalim, lolim);
    set(ax, 'Visible', 'off')
    states = shaperead('usastatehi', 'UseGeoCoords', true);
    geoshow(ax, states, 'FaceColor', [0.5 0.5 1])
	surfacem(xi,yi,aspect);
	%caxis([2 5]);
	colorbar;
	plotm(coor(:,1),coor(:,2),'rv');
	title('Surface Surface Gradient Direction');
	filename = sprintf('%s_%1d_azimap',event,period);
	print('-dpng',filename);
end
event=fgetl(fp);
end
