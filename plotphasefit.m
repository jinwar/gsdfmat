eventslist='eventlist';

lalim=[27 50.5];
lolim=[-114 -95];

gridsize=0.1;
lolim=lolim+360;
xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);

fp=fopen(eventslist,'r');

event=fgetl(fp);

while ischar(event)
for period=1:1

clear coor data ptxyz xyz xi yi z interpz;

% Read in the event la and lo
filename = sprintf('%s.log',event);
logfp=fopen(filename,'r');
stemp=fgetl(logfp);
fclose(logfp);
ftemp=sscanf(stemp,'Event: %f %f %f\n');
evla=ftemp(2);
evlo=ftemp(3);

% Read in the station arrival time
filename = sprintf('%s/best_%1d.stadt',event,period);

data= load(filename);

data(:,2)=data(:,2)+360;

% Kick out the station outside the boundary
[m n]=size(data)
j=1;
stalalim(1)=lalim(1)+2*gridsize;
stalalim(2)=lalim(2)-2*gridsize;
stalolim(1)=lolim(1)+2*gridsize;
stalolim(2)=lolim(2)-2*gridsize;
for i=1:m
	if data(i,1) < stalalim(2) && data(i,1) > stalalim(1) ...
			&& data(i,2) < stalolim(2) && data(i,2) > stalolim(1)
		coor(j,1)=data(i,1);
		coor(j,2)=data(i,2);
		coor(j,3)=data(i,3);
		coor(j,4)=distance(coor(j,1),coor(j,2),evla,evlo)*111/4;
		j=j+1;
	end
end

%ynode=ynode.*cosd(mean(xnode));
%coor(:,2)=coor(:,2).*cosd(mean(xnode));
[z,xi,yi]=gridfit(coor(:,1),coor(:,2),coor(:,3),xnode,ynode,...
					'smooth',2,'solver','normal');
%yi=yi./cosd(mean(xnode));
%coor(:,2)=coor(:,2)./cosd(mean(xnode));
%[testz,xi,yi]=gridfit(coor(:,1),coor(:,2),coor(:,4),xnode,ynode,...
%					'smooth',2,'solver','normal');

testz=griddata(coor(:,1),coor(:,2),coor(:,3),xi,yi);


figure(1);
clf
	usamap(lalim,lolim);
	[topoZ, refvec] = etopo('ETOPO5.DAT',1,lalim,lolim);
	geoshow(topoZ, refvec,'DisplayType','texturemap');
	colormap(demcmap(topoZ));
	plotm(coor(:,1),coor(:,2),'rv');
	title('Topography');
	filename = sprintf('%s_topo',event);
	print('-depsc',filename);

figure(2);
clf
    ax = usamap(lalim, lolim);
    set(ax, 'Visible', 'off')
    states = shaperead('usastatehi', 'UseGeoCoords', true);
    geoshow(ax, states, 'FaceColor', [0.5 0.5 1])
	contourm(xi,yi,z,30);
	colorbar;
	plotm(coor(:,1),coor(:,2),'rv');
    filename=sprintf('Phase Surface event:%s PID:%d',event,period);
	title(filename);
	filename = sprintf('%s_%1d_phmap',event,period);
	print('-depsc',filename);

% Calculate the gradient map
[aspect, slope, gradN, gradE]=gradientm(xi,yi,z);
[testaspect, testslope, testgradN, testgradE]=gradientm(xi,yi,testz);

% change the unit from slope to velocity
dz=tand(slope).^-1/1e3;
testdz=tand(testslope).^-1/1e3;

% Get rid of expolaration area
[m,n]=size(testdz);

% Calculate the azimuth of each grid
for i=1:m
	for j=1:n
		azi(i,j)=azimuth(xi(i,j), yi(i,j),evla,evlo);
	end
end
aspectch=azi;
for i=1:m
    for j=1:n
        %if abs(testdz(i,j)-4) < 0.05 && abs(testaspect(i,j)-azi(i,j)) < 1
		if testz(i,j)~=nan
			aspectch(i,j)=aspect(i,j)-azi(i,j);
		else
			aspectch(i,j)=nan;
			dz(i,j)=nan;
        end
    end
end

% Calculate the wave front direction
%[bigxi,bigyi]=meshgrid(lalim(1):(xnode(2)-xnode(1))*5:lalim(2), ...
%						lolim(1):(ynode(2)-ynode(1))*5:lolim(2)); 
%u=interp2(xi,yi,gradN,bigxi,bigyi);
%v=interp2(xi,yi,gradE,bigxi,bigyi);
%dz=tand(slope);

figure(3);
clf
    ax = usamap(lalim, lolim);
    set(ax, 'Visible', 'off')
    states = shaperead('usastatehi', 'UseGeoCoords', true);
    geoshow(ax, states, 'FaceColor', [0.5 0.5 1])
	surfacem(xi,yi,dz);
	seiscolormap
	load seiscmap
	colormap(seiscmap);
	%quiverm(bigxi,bigyi,u,v);
	%caxis([2 5]);
	colorbar;
	plotm(coor(:,1),coor(:,2),'rv');
    filename=sprintf('Surface Surface Gradient event:%s PID:%d',event,period);
    title(filename);
	filename = sprintf('%s_%1d_grdmap',event,period);
	print('-depsc',filename);

% Calculate the aspect change compare to the event back azimuth
figure(4);
clf
    ax = usamap(lalim, lolim);
    set(ax, 'Visible', 'off')
    states = shaperead('usastatehi', 'UseGeoCoords', true);
    geoshow(ax, states, 'FaceColor', [0.5 0.5 1])
	surfacem(xi,yi,aspectch);
	%caxis([2 5]);
	colorbar;
	plotm(coor(:,1),coor(:,2),'rv');
    filename=sprintf('Surface Surface Gradient Direction event:%s PID:%d',event,period);
    title(filename);
	filename = sprintf('%s_%1d_azimap',event,period);
	print('-depsc',filename);
end
event=fgetl(fp);
end
