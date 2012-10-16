% This script is to read in the output of program gsdfmain, calculate fit the 
% smooth phase surface, calculate the gradient, make the amplitude correction
% and output the final phase map for each event, each period.
%
% New in version 5: put in individual station amplitude check by comparing it 
% 	with surrounding stations.
%
% New in verson 4: put gmt surface fitting into consideration
%
% Ge Jin, jinwar@gmail.com
clear

% Input event list file
eventslist='testevent';

% some constants
ERRTOR=0.5;			 % the error allowed for cs measurement
ampmintor=0.3;
ampstator=0.5;
isonefigure=0;
phvrange=[[3 4];[3.5 4.5];[3.5 4.5]];
% phvrange(1,:)=[3.55 4.15];
periods=[25 50 100];


fp=fopen(eventslist,'r');

event=fgetl(fp);

while ischar(event)
for period=0:0

event
period
csinvfile=sprintf('%s_%1d.csinv',event,period);
stainvfile=sprintf('%s.stainv',event);
logfile=sprintf('%s.log',event);

lalim=[27 50.5];
lolim=[-115 -95];
lalim=[30 50];
lolim=[-125 -105];
% %lalim=[27 50.5];
%  lalim=[32 48];
%  lolim=[-113 -100];
gridsize=0.1;
lolim=lolim+360;
lat0=mean(lalim);
lon0=mean(lolim);
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
[stemp1 stemp2 evla evlo]=strread(stemp,'%s %s %f %f\n');

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

nerr=0;
for i=1:csnum
	err=csdata(i,3)-tnet(round(csdata(i,1))+1)+tnet(round(csdata(i,2)+1));
	if abs(err) > ERRTOR
		nerr=nerr+1;
	end
end
disp('Bad Measurement Number:');
disp(nerr);

% Iteratively exclude the bad measurements
oldcsdata=csdata;
newcsdata=csdata;
niter=1;
while nerr > 1 && niter<10
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
		if abs(err)> ERRTOR
			nerr=nerr+1;
		end
	end

	disp('Bad Measurement Number:');
	disp(nerr);

	oldcsdata=newcsdata;
	niter=niter+1;

end


% Find bad station and kick them out.
badsta=find(diag(A)<5);
for i=1:length(badsta)
	tnet(badsta(i))=nan;
end

% Find the amplitude abnormal stations and kick them out
for i=1:length(tnet)
	if ~isnan(tnet(i))
		csi=find(newcsdata(:,1)==i-1);
		stai1=newcsdata(csi,2)+1;
		csi=find(newcsdata(:,2)==i-1);
		stai2=newcsdata(csi,1)+1;
		stai=[stai1 ;stai2];
		ampmean=mean(sqrt(stadata(stai,4+period)));
		ampsta=sqrt(stadata(i,4+period));
		if ampsta < ampmean*ampstator || ampsta > ampmean/ampstator
			tnet(i)=nan;
		end
	end
end

% Check the list of bad stations and kick them out
if exist('badsta.lst')
	fp=fopen('badsta.lst','r');
	staname=fgetl(fp);
	while ischar(staname)
		stemp=sprintf('grep %s %s.sta | awk ''{print $2}'' > tempfile ',staname,event);
		system(stemp);
		staid=textread('tempfile','%n');
		if ~isempty(staid)
			tnet(staid+1)=nan;
		end
		staname=fgetl(fp);
	end
end

% Put the good stations' measurement in the data set.
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

disp('Number of Good stations');
disp(j);

% Transform the data from lat-lon domain into x-y domain
xycoor=coor;
[xycoor(:,1),xycoor(:,2)]=latlon2xy(coor(:,1),coor(:,2),lat0,lon0);
xyxnode=xnode-lat0;
xyynode=(ynode-lon0)*cosd(lalim(1));

% Fit the surface
%[z,xi,yi]=gridfit(coor(:,1),coor(:,2),coor(:,3),xnode,ynode,...
%					'smooth',2,'solver','normal');
[xy_z,xy_xi,xy_yi]=gridfit(xycoor(:,1),xycoor(:,2),xycoor(:,3),xyxnode,xyynode,...
					'smooth',2,'solver','normal');
% Fit the amplitude
[xy_ampmap,xy_xi,xy_yi]=gridfit(xycoor(:,1),xycoor(:,2),ap,xyxnode,xyynode,...
					'smooth',2,'solver','normal');

% Transform back to the lat-lon domain
[xi, yi]=meshgrid(xnode, ynode);
[m n]=size(xy_z);
k=0;
for i=1:m
	for j=1:n
		k=k+1;
		xy_data(k,1)=xy_xi(i,j);
		xy_data(k,2)=xy_yi(i,j);
		xy_data(k,3)=xy_z(i,j);
		xy_data(k,4)=xy_ampmap(i,j);
	end
end
lalo_data=xy_data;
[lalo_data(:,1),lalo_data(:,2)]=xy2latlon(xy_data(:,1),xy_data(:,2),lat0,lon0);
z=griddata(lalo_data(:,1),lalo_data(:,2),lalo_data(:,3),xi,yi,'cubic',{'QJ'});
ampmap=griddata(lalo_data(:,1),lalo_data(:,2),lalo_data(:,4),xi,yi,'cubic',{'QJ'});



% Fitting the phase and amplitude surface using GMT
ptxyz(:,1)=coor(:,2)-360;
ptxyz(:,2)=coor(:,1);
ptxyz(:,3)=coor(:,3);
save 'phasetemp.xyz' ptxyz -ASCII
system('csh phasemap.gmt');
xyz=load('gmtphsurf.xyz');
xyz(:,1)=xyz(:,1)+360;
gmtz=griddata(xyz(:,2),xyz(:,1),xyz(:,3),xi,yi,'cubic',{'QJ'}); % fit the phase map

ptxyz(:,3)=ap;
save 'phasetemp.xyz' ptxyz -ASCII
system('csh phasemap.gmt');
xyz=load('gmtphsurf.xyz');
xyz(:,1)=xyz(:,1)+360;
gmtamp=griddata(xyz(:,2),xyz(:,1),xyz(:,3),xi,yi,'cubic',{'QJ'}); % fit the amplitude map

% Calculate the correction term
dAmp=del2m(xi,yi,ampmap);
gmtdAmp=del2m(xi,yi,gmtamp);

% Calculate the correction term
amp_term=-dAmp./ampmap./(2*pi/periods(period+1)).^2;
gmtamp_term=-gmtdAmp./gmtamp./(2*pi/periods(period+1)).^2;


% Smooth it.
amp_term=smoothmap(xi,yi,amp_term);
gmtamp_term=smoothmap(xi,yi,gmtamp_term);

%[testz,xi,yi]=gridfit(coor(:,1),coor(:,2),coor(:,4),xnode,ynode,...
%					'smooth',2,'solver','normal');
%fs=fit([coor(:,1),coor(:,2)],coor(:,3),'cubicinterp');
%interpz=fs(xi,yi);
interpz=griddata(coor(:,1),coor(:,2),coor(:,3),xi,yi,'cubic',{'QJ'}); % interp z to find extended area


% Calculate the gradient and reliable area
[aspect, slope, gradN, gradE]=gradientm(xi,yi,z);
dz=tand(slope).^-1/1e3;
dz_correct=(dz.^-2+amp_term).^-.5;
correctV=-dz_correct+dz;

[aspect, slope, gradN, gradE]=gradientm(xi,yi,gmtz);	%use gmt result
gmtdz=tand(slope).^-1/1e3;
gmtdz_correct=(gmtdz.^-2+gmtamp_term).^-.5;	% use the gmt result

gmtcorrectV=-gmtdz_correct+gmtdz;

% Calculate the azimuth of each grid
[m,n]=size(xi);
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
		if ~isnan(interpz(i,j))
			aspectch(i,j)=aspect(i,j)-azi(i,j);
		else
			aspectch(i,j)=nan;
			dz(i,j)=nan;
			dz_correct(i,j)=nan;
			correctV(i,j)=nan;
			ampmap(i,j)=nan;
			gmtdz(i,j)=nan;
			gmtdz_correct(i,j)=nan;
			gmtcorrectV(i,j)=nan;

        end
    end
end

% Get rid of the small amplitude area and less constrained area
ampmean=nanmean(nanmean(ampmap));
[m n]=size(xi);

for i=1:m
	for j=1:n
		if abs(gmtz(i,j)-z(i,j))>1
			aspectch(i,j)=nan;
			dz(i,j)=nan;
			dz_correct(i,j)=nan;
			correctV(i,j)=nan;
			gmtdz(i,j)=nan;
			gmtdz_correct(i,j)=nan;
			gmtcorrectV(i,j)=nan;
			ampmap(i,j)=nan;
		end
	end
end

for i=1:m
	for j=1:n
		if ampmap(i,j)<ampmean*ampmintor
			aspectch(i,j)=nan;
			dz(i,j)=nan;
			dz_correct(i,j)=nan;
			correctV(i,j)=nan;
			gmtdz(i,j)=nan;
			gmtdz_correct(i,j)=nan;
			gmtcorrectV(i,j)=nan;
% 			ampmap(i,j)=nan;
		end
	end
end



%% Calculate the error of surface fitting
[m n]=size(xi);

	%find out the grids where the stations are
for k=1:length(coor)
	staj=find(abs(xi(1,:)-coor(k,1))<gridsize/2);
	if length(staj)>1
		staj=staj(1);
    end
	stai=find(abs(yi(:,1)-coor(k,2))<gridsize/2);
	if length(stai)>1
		stai=stai(1);
    end
    if isempty(stai) || isempty(staj)
        continue;
    end
	stamatrix(stai,staj)=1;
	staerr(k)=z(stai,staj)-coor(k,3);
	gmtstaerr(k)=gmtz(stai,staj)-coor(k,3);
end
%
%	% check the azi distribution of nearby station
%	% azicheck array is difined by:
%	% 2 | 1
%	% --+--
%	% 3 | 4
%	%
stad=2; % check range by degree
[m n]=size(xi);
relmatrix=zeros(m,n);
for i=1:m
	for j=1:n
		azicheck(1:4)=0;
		if ~isnan(dz(i,j))
			staj=find(abs(coor(:,1)-xi(i,j))<stad);
			stai=find(abs(coor(staj,2)-yi(i,j))<stad/cosd(xi(i,j)));
			staj=staj(stai);
			stai=find(dist(coor(staj,1),coor(staj,2),xi(i,j),yi(i,j))<stad*111);
			stai=staj(stai);
			if ~isempty(stai)
				for k=1:length(stai)
					dla=coor(stai(k),1)-xi(i,j);
					dlo=coor(stai(k),2)-yi(i,j);
					dlo=dlo*cosd(xi(i,j));
					if dla>0 && dla<stad
						if dlo>0 && dlo<stad
							azicheck(1)=1;
						elseif dlo>-stad
							azicheck(2)=1;
						end
					elseif dla>-stad
						if dlo>0 && dlo<stad
							azicheck(4)=1;
						elseif dlo>-stad
							azicheck(3)=1;
						end
					end
				end
%				if sum(azicheck(1:4))>3.9
%					relmatrix(i,j)=1;
%				end
				relmatrix(i,j)=sum(azicheck(1:4));
			end
		end
	end
end


for i=1:m
	for j=1:n
		if relmatrix(i,j)<3.5 && ~isnan(dz(i,j))
			aspectch(i,j)=nan;
			dz(i,j)=nan;
			dz_correct(i,j)=nan;
			correctV(i,j)=nan;
			ampmap(i,j)=nan;
		end
	end
end

for i=1:m
    for j=1:n
        if isnan(dz(i,j))
            z(i,j)=nan;
        end
    end
end

% save the workspace
filename = sprintf('%s_%1d.mat',event,period);
save(filename)


if isonefigure
	drawevent(event,period);
	figure(2)
	clf
		ax = usamap(lalim, lolim);
		set(ax, 'Visible', 'off')
		states = shaperead('usastatehi', 'UseGeoCoords', true);
		geoshow(ax, states, 'FaceColor', [0.5 0.5 1])
		surfacem(xi,yi,gmtdz_correct);
		contourm(xi,yi,gmtz,30,'k');
		geoshow(ax, states, 'FaceColor', 'none')
		seiscolormap
		load seiscmap
		colormap(seiscmap);
		%quiverm(bigxi,bigyi,u,v);
		caxis(phvrange(period+1,:));
		colorbar;
		plotm(coor(:,1),coor(:,2),'rv');
else

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
		plotm(coor(:,1),coor(:,2),'rv');
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
        geoshow(ax, states, 'FaceColor', [0.5 0.5 1],'linewidth',1.5)
		surfacem(xi,yi,ampmap);
        plotm(coor(:,1),coor(:,2),'rv');
		%caxis([2 5]);
		colorbar;
        set(gca,'fontsize',15)
% 		plotm(coor(:,1),coor(:,2),'rv');
		filename=sprintf('Amplitude Map event:%s PID:%d',event,period);
% 		title(filename,'Interpreter','none');
		filename = sprintf('%s_%1d_ampmap',event,period)
		print('-depsc',filename);

	figure(4)
	clf
		ax = usamap(lalim, lolim);
		set(ax, 'Visible', 'off')
		states = shaperead('usastatehi', 'UseGeoCoords', true);
        geoshow(ax, states, 'FaceColor', [0.5 0.5 1],'linewidth',1.5)
        surfacem(xi,yi,z);
				contourm(xi,yi,z,10,'k','linewidth',1.5);
%         sc=nanmean(nanstd(correctV));
% 		caxis([-4*sc 4*sc]);
% 		seiscolormap
% 		load seiscmap
% 		colormap(seiscmap);
		colorbar;
		plotm(coor(:,1),coor(:,2),'rv');
		filename=sprintf('Amplitude Correction Map event:%s PID:%d',event,period);
% 		title(filename,'Interpreter','none');
% 		filename = sprintf('%s_%1d_ampcorrection',event,period)
 		filename = sprintf('%s_%1d_phasemap',event,period)
		print('-depsc',filename);

	figure(5);
	clf
%         lalim=[30 48];
%         lolim=[-113 -98];
		ax = usamap(lalim, lolim);
		set(ax, 'Visible', 'off')
		states = shaperead('usastatehi', 'UseGeoCoords', true);
		geoshow(ax, states, 'FaceColor', [0.5 0.5 1],'linewidth',1.5)
		surfacem(xi,yi,dz_correct);
		contourm(xi,yi,z,10,'k','linewidth',1.5);
		geoshow(ax, states, 'FaceColor', 'none')
		seiscolormap
		load seiscmap
		colormap(seiscmap);
		%quiverm(bigxi,bigyi,u,v);
		caxis(phvrange(period+1,:));
		colorbar;
        set(gca,'fontsize',15)
% 		plotm(coor(:,1),coor(:,2),'rv');
		filename = sprintf('%s_%1d_grdmap_correct',event,period)
% 		title(filename,'Interpreter','none');
		print('-depsc',filename);

end % End of figure plot if


end	% end of period loop
event=fgetl(fp);

end % end of event loop
