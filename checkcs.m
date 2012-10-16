clear;

system('csh makesimcs.csh');
data=load('200910020107_2.simcs');
evla=-16.43;
evlo=-173.05;

lalim=[27 50.5];
lolim=[-114 -95];


[dist1, azi1]=distance(evla,evlo,data(:,1),data(:,2));
dist1=deg2km(dist1);
for i=1:length(dist1)
	dist1(i)=vdist(evla, evlo, data(i,1),data(i,2))/1e3;
end
[dist2, azi2]=distance(evla,evlo,data(:,3),data(:,4));
dist2=deg2km(dist2);
for i=1:length(dist1)
	dist2(i)=vdist(evla, evlo, data(i,3),data(i,4))/1e3;
end

ddist=dist1-dist2;
dazi=azi1-azi2;
[sdist sazi]=distance(data(:,1),data(:,2),data(:,3),data(:,4));
sdist=deg2km(sdist);
avgdist=(dist1+dist2)/2;
dt=data(:,5);

v=polyfit(dt,ddist,1)
dterr=dt-ddist/v(1);

badindex=find(abs(dterr)>0.2);

for i=1:length(badindex)
	stemp=sprintf('grep %6.3f 200910020107.cs | grep %6.3f |awk ''{if ($10 == 100 && $16 < 5) print $1" "$2" "%f}'' ',data(badindex(i),1),data(badindex(i),3),dterr(badindex(i)));
	%system(stemp);
end


figure(1)
clf;
ax = usamap(lalim, lolim);
set(ax, 'Visible', 'off')
states = shaperead('usastatehi', 'UseGeoCoords', true);
geoshow(ax, states, 'FaceColor', [0.5 0.5 1])
hold on
for i=1:length(dterr)
	stala=[data(i,1),data(i,3)];
	stalo=[data(i,2),data(i,4)]+360;
	if abs(dterr(i))>0.2
		s='r';
	elseif abs(dterr(i))>0.1
		s='b';
	else
		s='.';
	end
	plotm(stala,stalo,s);
end

for i=find(abs(dterr)>0.1)
	stala=[data(i,1),data(i,3)];
	stalo=[data(i,2),data(i,4)]+360;
	if abs(dterr(i))>0.2
		s='r';
	elseif abs(dterr(i))>0.1
		s='b';
	else
		s='.';
	end
	plotm(stala,stalo,s);
end

figure(3)
clf;
plot(ddist,dterr,'.')

figure(2)
clf
plot(dazi,dterr,'.');

figure(4)
clf
plot(sazi,dterr,'.');

