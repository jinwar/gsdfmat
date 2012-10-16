% This script is to read in the output of program gsdfmain, calculate fit the 
% smooth phase surface, calculate the 
clear

eventslist='goodevent';
fp=fopen(eventslist,'r');

event=fgetl(fp);

while ischar(event)
for period=0:2

csinvfile=sprintf('%s_%1d.csinv',event,period);
stainvfile=sprintf('%s.stainv',event);
logfile=sprintf('%s.log',event);
phvrange=[[3 4];[3.5 4.5];[3.5 4.5]];
periods=[25 50 100];

lalim=[27 50.5];
lolim=[-114 -95];
gridsize=0.1;
lolim=lolim+360;
xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);

% read in data and information
csdata=load(csinvfile);
stadata=load(stainvfile);
stanum=length(stadata);
csnum=length(csdata);

logfp=fopen(logfile,'r');
stemp=fgetl(logfp);
fclose(logfp);
ftemp=sscanf(stemp,'Event: %f %f %f\n');
evla=ftemp(2);
evlo=ftemp(3);

% Make the matrix
mat=sparse(csnum,stanum);
dt=zeros(csnum,1);

for i=1:csnum
	mat(i,round(csdata(i,1))+1)=1;
	mat(i,round(csdata(i,2))+1)=-1;
	dt(i)=csdata(i,3);
end

% find the best station as reference station
A=mat'*mat;
maxdiag=max(diag(A));
beststa=find(diag(A)==maxdiag);
mat(csnum+1,beststa(1))=1;
dt(csnum+1)=0;

A=mat'*mat;
disp('start inverse');
tnet=A\(mat'*dt);

% Exclude bad measurement
ERRTOR=0.5;

nerr=0;
for i=1:csnum
	err=csdata(i,3)-tnet(round(csdata(i,1))+1)+tnet(round(csdata(i,2)+1));
	if abs(err) > ERRTOR
		nerr=nerr+1;
	end
end
disp('Bad Measurement Number:');
disp(nerr);

oldcsdata=csdata;
newcsdata=csdata;
while nerr > 10
	clear mat dt newcsdata;
	mat=sparse(length(oldcsdata)-nerr,stanum);
	dt=zeros(length(oldcsdata)-nerr,1);

	j=1;
	for i=1:length(oldcsdata)
		err=oldcsdata(i,3)-tnet(round(oldcsdata(i,1))+1)+tnet(round(oldcsdata(i,2)+1));
		if abs(err) <= ERRTOR
			newcsdata(j,:)=oldcsdata(i,:);
			mat(j,round(oldcsdata(i,1))+1)=1;
			mat(j,round(oldcsdata(i,2))+1)=-1;
			dt(j)=oldcsdata(i,3);
            j=j+1;
		end
	end

	% find the best station as reference station
	maxdiag=max(diag(A));
	beststa=find(diag(A)==maxdiag);
	mat(csnum+1,beststa(1))=1;
	dt(csnum+1)=0;

	A=mat'*mat;
	disp('start inverse');
	tnet=A\(mat'*dt);

	nerr=0;
	for i=1:length(newcsdata)
		err=newcsdata(i,3)-tnet(round(newcsdata(i,1))+1)+tnet(round(newcsdata(i,2)+1));
		if abs(err)>1
			nerr=nerr+1;
		end
	end

	disp('Bad Measurement Number:');
	disp(nerr);

	oldcsdata=newcsdata;

end


% Find bad station and kick them out.
badsta=find(diag(A)<5);
for i=1:length(badsta)
	tnet(badsta(i))=nan;
end

coor=0;ap=0;
clear coor,ap;
j=1;
for i=1:length(tnet)
	if ~isnan(tnet(i))
		coor(j,1:2)=stadata(i,2:3);
		coor(j,2)=coor(j,2)+360;
		coor(j,3)=tnet(i);
		coor(j,4)=distance(coor(j,1),coor(j,2),evla,evlo)*111/4;
		ap(j)=sqrt(stadata(i,4+period));
		j=j+1;
	end
end

% Fit the surface
[z,xi,yi]=gridfit(coor(:,1),coor(:,2),coor(:,3),xnode,ynode,...
					'smooth',2,'solver','normal');
% Fit the amplitude
[ampmap,xi,yi]=gridfit(coor(:,1),coor(:,2),ap,xnode,ynode,...
					'smooth',2,'solver','normal');

% Calculate the correction term
dAmp=del2m(xi,yi,ampmap);

% Calculate the correction term
amp_term=-dAmp./ampmap./(2*pi/periods(period+1)).^2;

% Smooth it.
amp_term=smoothmap(xi,yi,amp_term);

%[testz,xi,yi]=gridfit(coor(:,1),coor(:,2),coor(:,4),xnode,ynode,...
%					'smooth',2,'solver','normal');
testz=griddata(coor(:,1),coor(:,2),coor(:,3),xi,yi);


% Calculate the gradient and reliable area
[aspect, slope, gradN, gradE]=gradientm(xi,yi,z);
[testaspect, testslope, testgradN, testgradE]=gradientm(xi,yi,testz);
dz=tand(slope).^-1/1e3;
testdz=tand(testslope).^-1/1e3;
dz_correct=(dz.^-2+amp_term).^-.5;
correctV=-dz_correct+dz;

% Calculate the azimuth of each grid
[m,n]=size(testdz);
for i=1:m
	for j=1:n
		azi(i,j)=azimuth(xi(i,j), yi(i,j),evla,evlo);
	end
end

% Get rid of expolaration area
aspectch=azi;
for i=1:m
    for j=1:n
        %if abs(testdz(i,j)-4) < 0.05 && abs(testaspect(i,j)-azi(i,j)) < 1
		if ~isnan(testz(i,j))
			aspectch(i,j)=aspect(i,j)-azi(i,j);
		else
			aspectch(i,j)=nan;
			dz(i,j)=nan;
			dz_correct(i,j)=nan;
			correctV(i,j)=nan;
        end
    end
end

figure(1);
clf
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
	print('-depsc',filename);

% Calculate the aspect change compare to the event back azimuth
figure(2);
clf
    ax = usamap(lalim, lolim);
    set(ax, 'Visible', 'off')
    states = shaperead('usastatehi', 'UseGeoCoords', true);
    geoshow(ax, states, 'FaceColor', [0.5 0.5 1])
	surfacem(xi,yi,aspectch);
	%contourm(xi,yi,z,30,'k');
	geoshow(ax, states, 'FaceColor', 'none')
	%caxis([2 5]);
	colorbar;
	plotm(coor(:,1),coor(:,2),'rv');
    filename=sprintf('Surface Surface Gradient Direction event:%s PID:%d',event,period);
	title(filename,'Interpreter','none');
	filename = sprintf('%s_%1d_azimap',event,period);
	print('-depsc',filename);

figure(3);
clf
    ax = usamap(lalim, lolim);
    set(ax, 'Visible', 'off')
    states = shaperead('usastatehi', 'UseGeoCoords', true);
    geoshow(ax, states, 'FaceColor', [0.5 0.5 1])
	surfacem(xi,yi,ampmap);
	%caxis([2 5]);
	colorbar;
	plotm(coor(:,1),coor(:,2),'rv');
    filename=sprintf('Amplitude Map event:%s PID:%d',event,period);
	title(filename,'Interpreter','none');
	filename = sprintf('%s_%1d_ampmap',event,period);
	print('-depsc',filename);

figure(4)
clf
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
	print('-depsc',filename);

figure(5);
clf
	ax = usamap(lalim, lolim);
	set(ax, 'Visible', 'off')
	states = shaperead('usastatehi', 'UseGeoCoords', true);
	geoshow(ax, states, 'FaceColor', [0.5 0.5 1])
	surfacem(xi,yi,dz_correct);
	contourm(xi,yi,z,30,'k');
	geoshow(ax, states, 'FaceColor', 'none')
	seiscolormap
	load seiscmap
	colormap(seiscmap);
	%quiverm(bigxi,bigyi,u,v);
	caxis(phvrange(period+1,:));
	colorbar;
	%plotm(coor(:,1),coor(:,2),'rv');
	filename = sprintf('%s_%1d_grdmap_correct',event,period);
	title(filename,'Interpreter','none');
	print('-depsc',filename);


% save the output
filename = sprintf('%s_%1d.mat',event,period);
save(filename)

end	% end of period loop
event=fgetl(fp);

end % end of event loop
