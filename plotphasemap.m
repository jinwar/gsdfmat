event='200601040832'
period=1;

filename = sprintf('%s/best_%1d.stadt',event,period);

data= load(filename);

coor=data(:,1:2);
coor(:,2)=coor(:,2)+360;

lalim=[32 43.5];
lolim=[-125 -112];
lolim=lolim+360;

[xi,yi]=meshgrid(lalim(1):0.1:lalim(2), lolim(1):0.1:lolim(2)); 
z=griddata(coor(:,1),coor(:,2),data(:,3),xi,yi,'cubic');

[x,y,z]=grdread2('image.bin');

[xi,yi]=meshgrid(lalim(1):0.1:lalim(2), lolim(1):0.1:lolim(2));

figure(1);
clf
	usamap(lalim,lolim);
	[Z, refvec] = etopo('ETOPO5.DAT',1,lalim,lolim);
	geoshow(Z, refvec,'DisplayType','texturemap');
	colormap(demcmap(Z));
	plotm(coor(:,1),coor(:,2),'rv');

figure(2);
clf
    ax = usamap(lalim, lolim);
    set(ax, 'Visible', 'off')
    states = shaperead('usastatehi', 'UseGeoCoords', true);
    geoshow(ax, states, 'FaceColor', [0.5 0.5 1])
	surfacem(xi,yi,z);
	colorbar;
	plotm(coor(:,1),coor(:,2),'rv');

[aspect, slope, gradN, gradE]=gradientm(xi,yi,z);
dz=tand(slope).^-1/1e3;
%dz=tand(slope);

figure(3);
clf
    ax = usamap(lalim, lolim);
    set(ax, 'Visible', 'off')
    states = shaperead('usastatehi', 'UseGeoCoords', true);
    geoshow(ax, states, 'FaceColor', [0.5 0.5 1])
	surfacem(xi,yi,dz);
	%caxis([2 5]);
	colorbar;
	plotm(coor(:,1),coor(:,2),'rv');

