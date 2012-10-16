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
filename = sprintf('%s/N26A_%1d.stadt',event,period);

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
[testz,xi,yi]=gridfit(coor(:,1),coor(:,2),coor(:,4),xnode,ynode,...
					'smooth',2,'solver','normal');


figure(1);
clf
	usamap(lalim,lolim);
	[topoZ, refvec] = etopo('ETOPO5.DAT',1,lalim,lolim);
	geoshow(topoZ, refvec,'DisplayType','texturemap');
	colormap(demcmap(topoZ));
	plotm(coor(:,1),coor(:,2),'rv');
	title('Topography');
	filename = sprintf('%s_topo',event);
	print('-dpng',filename);

figure(2);
clf
    ax = usamap(lalim, lolim);
    set(ax, 'Visible', 'off')
    states = shaperead('usastatehi', 'UseGeoCoords', true);
	geoshow(xi,yi,z,'DisplayType','texture');
	contourm(xi,yi,z,30,'b');
    geoshow(ax, states,'facecolor','none')
	colorbar;
	plotm(coor(:,1),coor(:,2),'kv');
    filename=sprintf('Event:%s Frequency Band:50s',event);
	title(filename);
	filename = sprintf('%s_%1d_phmap',event,period);
	print('-djpeg99',filename);

end
event=fgetl(fp);
end
