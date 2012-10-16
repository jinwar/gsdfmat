% This script is to read in the output of program gsdfmain, calculate fit the 
% smooth phase surface, calculate the gradient, make the amplitude correction
% and output the final phase map for each event, each period.
%
% New in the version 7: change the method to build up the smoothing kernel, from gaussian smoothing to Laplician Term
%
% New in version 6: instead of using surface fitting to calculate the phase surface
% then make the gradient, we inverse for the phase gradient surface
%
% New in version 5: put in individual station amplitude check by comparing it 
% 	with surrounding stations.
%
% New in verson 4: put gmt surface fitting into consideration
%
% Ge Jin, jinwar@gmail.com
clear

% Input event list file
eventslist='allcs';

% some constants
ERRTOR=0.5;			 % the error allowed for cs measurement
mincsnum=100;
isfigure=0;

ampmintor=0.01;
ampstator=0.1;
isfigure=0;
avgphv=[3.53 3.64 3.73 3.78 3.83 3.94 4.06 4.27];
r=0.05;
phvrange(:,1)=avgphv'.*(1-r);
phvrange(:,2)=avgphv'.*(1+r);

%phvrange(1,:)=[3.55 4.15];
periods=[20 25 32 40 50 66 83 100];

periodrange=0:7;
smweight0 = 2;


lalim=[-11.2 -7.8];
lolim=[148.8 151.5];
gridsize=0.1;

%lalim=[30 50];
%lolim=[-125 -90];
%gridsize=0.3;

raydensetol=deg2km(gridsize)*1;
lolim=lolim;
lat0=mean(lalim);
lon0=mean(lolim);
xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);

Nx=length(xnode);
Ny=length(ynode);
[xi yi]=ndgrid(xnode,ynode);

% new, using Laplacian term
disp('initial the smoothing kernel')
tic
    [i,j] = ndgrid(1:Nx,2:(Ny-1));
    ind = j(:) + Ny*(i(:)-1);
    dy = diff(ynode);
    dy1 = dy(j(:)-1);
    dy2 = dy(j(:));

    Areg = sparse(repmat(ind,1,3),[ind-1,ind,ind+1], ...
                    [-2./(dy1.*(dy1+dy2)), 2./(dy1.*dy2), -2./(dy2.*(dy1+dy2))],Nx*Ny,Nx*Ny);

    [i,j] = ndgrid(2:(Nx-1),1:Ny);
    ind = j(:) + Ny*(i(:)-1);
    dx = diff(xnode);
    dx1 = dx(i(:)-1);
    dx2 = dx(i(:));

    Areg = [Areg;sparse(repmat(ind,1,3),[ind-Ny,ind,ind+Ny], ...
            [-2./(dx1.*(dx1+dx2)), 2./(dx1.*dx2), -2./(dx2.*(dx1+dx2))],Nx*Ny,Nx*Ny)];

    F=sparse(Nx*Ny*2*2,Nx*Ny*2);
    oldprocess=0;
    for n=1:size(Areg,1)
	   newprocess=floor(n/2/Nx/Ny*10);
	   if newprocess>oldprocess
		   disp(['Finished: ', num2str(newprocess), '0%'])
		   oldprocess=newprocess;
	   end
        ind=find(Areg(n,:)~=0);
        F(2*n-1,2*ind-1)=Areg(n,ind);
        F(2*n,2*ind)=Areg(n,ind);
    end
toc


fp=fopen(eventslist,'r');

event=fgetl(fp);

while ischar(event)
for period=periodrange

event
period
csinvfile=sprintf('%s_%1d.csinv',event,period);
stainvfile=sprintf('%s.stainv',event);
stafile=sprintf('%s.sta',event);
logfile=sprintf('%s.log',event);


%lalim=[27 50.5];


% read in data and information
csdata=load(csinvfile);
stadata=load(stainvfile);
stanum=length(stadata);
csnum=size(csdata,1);
if csnum < mincsnum
	continue;
end

logfp=fopen(logfile,'r');
stemp=fgetl(logfp);
fclose(logfp);
[stemp1 stemp2 evla evlo]=strread(stemp,'%s %s %f %f\n');

% Make the matrix
%mat=sparse(csnum,stanum);
dt=zeros(csnum,1);
rays=zeros(csnum,4);
W = sparse(length(dt),length(dt));
for i=1:length(dt)
	W(i,i)=1;
end
for i=1:csnum
%	mat(i,round(csdata(i,1))+1)=1;
%	mat(i,round(csdata(i,2))+1)=-1;
    sta1coor=stadata(round(csdata(i,1))+1,2:3);
    sta2coor=stadata(round(csdata(i,2))+1,2:3);
    Isinmap=1;
    % Check whether it's in the range
    temp=[sta1coor(1) sta2coor(1)];
    if min(temp) < lalim(1) || max(temp) > lalim(2)
        Isinmap=0;
    end

    temp=[sta1coor(2) sta2coor(2)];
    if min(temp) < lolim(1) || max(temp) > lolim(2)
        Isinmap=0;
    end
    if Isinmap
        rays(i,1:2)=stadata(round(csdata(i,1))+1,2:3);
        rays(i,3:4)=stadata(round(csdata(i,2))+1,2:3);
        dt(i)=csdata(i,3);
    else
        rays(i,1:2)=0;
        rays(i,3:4)=0;
        dt(i)=0;
        W(i,i)=0;
    end
end
disp('Start building the kernel');
tic
mat=kernel_build(rays,xnode,ynode);
toc

% Calculate the kernel density
%sumG=sum(abs(mat),1);
ind=1:Nx*Ny;
sumG(ind)=sum((mat(:,2*ind).^2+mat(:,2*ind-1).^2).^.5,1);
for i=1:Nx
	for j=1:Ny
		n=Ny*(i-1)+j;
		raydense(i,j)=sumG(n);
	end
end

% build the smoothing operator

smweight = smweight0;
NR=norm(F,1);
NA=norm(mat,1);
smweight = smweight0*NA/NR;

disp('start inverse');

    A=[W*mat;smweight*F];
    rhs=[W*dtp;zeros(size(F,1),1)];

    disp('start inverse');
    tic
		phaseg=(A'*A)\(A'*rhs);
    toc
    disp('Done');


% Iteratively down weight the measurement with high error
niter=1;

while niter < 2
	niter
	niter=niter+1;
	err = mat*phaseg - dt;
	err = W*err;
	stderr=std(err);
	for i=1:length(err)
		if abs(err(i)) > 2*stderr
			W(i,i)=0;
		end
	end

	% Rescale the smooth kernel
    NR=norm(F,1);
    NA=norm(mat,1);
    smweight = smweight0*NA/NR;

	A=[W*mat;smweight*F];
	rhs=[W*dtp;zeros(size(F,1),1)];

	disp('start inverse');
	tic
		phaseg=(A'*A)\(A'*rhs);
	toc
	disp('Done');

end

disp(' Get rid of uncertainty area');
%sumG=sum(abs(mat),1);
for i=1:Nx
	for j=1:Ny
		n=Ny*(i-1)+j;
%		raydense(i,j)=(sumG(2*n-1).^2 + sumG(2*n).^2).^5;
		if raydense(i,j) < raydensetol
			phaseg(2*n-1)=NaN;
			phaseg(2*n)=NaN;
		end
	end
end

%plotphaseg(phaseg,stadata,xnode,ynode);

% Change phaseg into phase velocity
for i=1:Nx
	for j=1:Ny
		n=Ny*(i-1)+j;
		GVx(i,j)= phaseg(2*n-1);
		GVy(i,j)= phaseg(2*n);
	end
end

GV=(GVx.^2+GVy.^2).^-.5;

% Find the amplitude abnormal stations and kick them out
goodcsnum=0;
clear newcsdata;
for i=1:size(csdata,1)
	if W(i,i)==1
		goodcsnum=goodcsnum+1;
		newcsdata(goodcsnum,:)=csdata(i,:);
	end
end

for i=1:size(stadata,1)
    if exist('newcsdata','var')
        csi=find(newcsdata(:,1)==i-1);
        stai1=newcsdata(csi,2)+1;
        csi=find(newcsdata(:,2)==i-1);
        stai2=newcsdata(csi,1)+1;
        stai=[stai1 ;stai2];
    else
        stai=[];
    end
	if isempty(stai)
		IsGoodSta(i)=0;
	else
		ampmean=mean(sqrt(stadata(stai,4+period)));
		ampsta=sqrt(stadata(i,4+period));
		if ampsta < ampmean*ampstator || ampsta > ampmean/ampstator
			IsGoodSta(i)=0;
		else
			IsGoodSta(i)=1;
		end
    end
end

for i=1:size(stadata,1)
    if stadata(i,2)<lalim(1) || stadata(i,2)>lalim(2)
        IsGoodSta(i)=0;
    end
    if stadata(i,3)<lolim(1) || stadata(i,3)>lolim(2)
        IsGoodSta(i)=0;
    end
end

% Check the list of bad stations and kick them out
if exist('badsta.lst')
	fp_badsta=fopen('badsta.lst','r');
	staname=fgetl(fp_badsta);
	while ischar(staname)
		stemp=sprintf('grep %s %s.sta | awk ''{print $2}'' > tempfile ',staname,event);
		system(stemp);
		staid=textread('tempfile','%n');
		if ~isempty(staid)
			IsGoodSta(staid+1)=0;
		end
		staname=fgetl(fp_badsta);
	end
	fclose(fp_badsta);
end

% Put the good stations' measurement in the data set.
coor=0;ap=0;
clear coor,ap;
j=1;
for i=1:size(stadata,1)
	if IsGoodSta(i)
		coor(j,1:2)=stadata(i,2:3);
		coor(j,2)=coor(j,2);
		coor(j,3)=distance(coor(j,1),coor(j,2),evla,evlo)*111/4;
		ap(j)=sqrt(stadata(i,4+period));
		j=j+1;
	else
% 		fp_sta=fopen(stafile,'r');
% 		for k=1:i+1
% 			stemp=fgetl(fp_sta);
% 		end
% 		disp(['Bad Station: ', stemp])
% 		fclose(fp_sta);
	end
end

disp(['Number of Bad stations:' num2str(size(stadata,1)-j)]);

disp(['Number of Good stations:' num2str(j) ]);

% Transform the data from lat-lon domain into x-y domain
xycoor=coor;
[xycoor(:,1),xycoor(:,2)]=latlon2xy(coor(:,1),coor(:,2),lat0,lon0);
xyxnode=xnode-lat0;
xyynode=(ynode-lon0)*cosd(lalim(1));

% Fit the surface
%[z,xi,yi]=gridfit(coor(:,1),coor(:,2),coor(:,3),xnode,ynode,...
%					'smooth',2,'solver','normal');
%[xy_z,xy_xi,xy_yi]=gridfit(xycoor(:,1),xycoor(:,2),xycoor(:,3),xyxnode,xyynode,...
%					'smooth',2,'solver','normal');
% Fit the amplitude
[xy_ampmap,xy_xi,xy_yi]=gridfit_jg(xycoor(:,1),xycoor(:,2),ap,xyxnode,xyynode,...
					'smooth',2,'regularizer','del4','solver','normal');

% Transform back to the lat-lon domain
[xi, yi]=meshgrid(xnode, ynode);
[m n]=size(xy_xi);
k=0;
for i=1:m
	for j=1:n
		k=k+1;
		xy_data(k,1)=xy_xi(i,j);
		xy_data(k,2)=xy_yi(i,j);
%		xy_data(k,3)=xy_z(i,j);
		xy_data(k,3)=xy_ampmap(i,j);
	end
end
lalo_data=xy_data;
[lalo_data(:,1),lalo_data(:,2)]=xy2latlon(xy_data(:,1),xy_data(:,2),lat0,lon0);
ampmap=griddata(lalo_data(:,1),lalo_data(:,2),lalo_data(:,3),xi,yi,'cubic',{'QJ'});

dAmp=del2m(xi,yi,ampmap);

% Calculate the correction term
amp_term=-dAmp./ampmap./(2*pi/periods(period+1)).^2;


% Smooth it.
smD=max([300 periods(period+1).*mean(phvrange(period+1,:))]);
amp_term=smoothmap(xi,yi,amp_term,smD);

% Correct the phase anormaly
GV_correct=((GV).^-2 + amp_term').^-.5;

if isfigure
	figure(1)
	clf
	hold on
		ax = worldmap(lalim, lolim);
		set(ax, 'Visible', 'off')
		states = shaperead('usastatehi', 'UseGeoCoords', true);
		surfacem(xi,yi,GV')
%		plotm(stadata(:,2),stadata(:,3),'rv');
		%contourm(xi,yi,z,30,'k');
		geoshow(ax, states, 'FaceColor', 'none')
		seiscolormap
		load seiscmap
		colormap(seiscmap);
		colorbar
		caxis(phvrange(period+1,:))
		filename = sprintf('%s_%1d_grdmap',event,period);
		title(filename,'Interpreter','none');
		print('-depsc',filename);

	figure(2)
	clf
	hold on
		ax = usamap(lalim, lolim);
		set(ax, 'Visible', 'off')
		states = shaperead('usastatehi', 'UseGeoCoords', true);
		surfacem(xi,yi,GV_correct')
%		plotm(stadata(:,2),stadata(:,3),'rv');
		%contourm(xi,yi,z,30,'k');
		geoshow(ax, states, 'FaceColor', 'none')
		seiscolormap
		load seiscmap
		colormap(seiscmap);
		colorbar
		caxis(phvrange(period+1,:))
		filename = sprintf('%s_%1d_grdmap_correct',event,period);
		title(filename,'Interpreter','none');
		print('-depsc',filename);

    figure(3)
	clf
	hold on
		ax = usamap(lalim, lolim);
		set(ax, 'Visible', 'off')
		states = shaperead('usastatehi', 'UseGeoCoords', true);
		surfacem(xi,yi,GV_correct'-GV')
%		plotm(stadata(:,2),stadata(:,3),'rv');
		%contourm(xi,yi,z,30,'k');
		geoshow(ax, states, 'FaceColor', 'none')
		seiscolormap
		load seiscmap
		colormap(seiscmap);
		colorbar
%		caxis([3.5 4.15])
		filename = sprintf('%s_%1d_correctterm',event,period);
		title(filename,'Interpreter','none');
		print('-depsc',filename);

	figure(4)
	clf
	hold on
		ax = usamap(lalim, lolim);
		set(ax, 'Visible', 'off')
		states = shaperead('usastatehi', 'UseGeoCoords', true);
		surfacem(xi,yi,ampmap)
		plotm(stadata(:,2),stadata(:,3),'rv');
		%contourm(xi,yi,z,30,'k');
		geoshow(ax, states, 'FaceColor', 'none')
%		seiscolormap
%		load seiscmap
%		colormap(seiscmap);
		colorbar
		filename = sprintf('%s_%1d_ampmap',event,period);
		title(filename,'Interpreter','none');
		print('-depsc',filename);
end
%[testz,xi,yi]=gridfit(coor(:,1),coor(:,2),coor(:,4),xnode,ynode,...
%					'smooth',2,'solver','normal');
%fs=fit([coor(:,1),coor(:,2)],coor(:,3),'cubicinterp');
%interpz=fs(xi,yi);
%interpz=griddata(coor(:,1),coor(:,2),coor(:,3),xi,yi,'cubic',{'QJ'}); % interp z to find extended area
%
%
%% Calculate the gradient and reliable area
%[aspect, slope, gradN, gradE]=gradientm(xi,yi,z);
%dz=tand(slope).^-1/1e3;
%correctV=-dz_correct+dz;

%[aspect, slope, gradN, gradE]=gradientm(xi,yi,gmtz);	%use gmt result
%gmtdz=tand(slope).^-1/1e3;
%gmtdz_correct=(gmtdz.^-2+gmtamp_term).^-.5;	% use the gmt result
%
%gmtcorrectV=-gmtdz_correct+gmtdz;
%
%% Calculate the azimuth of each grid
%[m,n]=size(xi);
%for i=1:m
%	for j=1:n
%		azi(i,j)=azimuth(xi(i,j), yi(i,j),evla,evlo);
%	end
%end
%
% Get rid of expolaration area
%aspectch=azi;
%for i=1:m
%    for j=1:n
%        %if abs(testdz(i,j)-4) < 0.05 && abs(testaspect(i,j)-azi(i,j)) < 1
%		if ~isnan(interpz(i,j))
%			aspectch(i,j)=aspect(i,j)-azi(i,j);
%		else
%			aspectch(i,j)=nan;
%			dz(i,j)=nan;
%			dz_correct(i,j)=nan;
%			correctV(i,j)=nan;
%			ampmap(i,j)=nan;
%			gmtdz(i,j)=nan;
%			gmtdz_correct(i,j)=nan;
%			gmtcorrectV(i,j)=nan;
%
%        end
%    end
%end

% Get rid of the small amplitude area 
ampmean=nanmean(nanmean(ampmap));
[m n]=size(xi);

for i=1:m
	for j=1:n
		if ampmap(i,j)<ampmean*ampmintor
%			aspectch(i,j)=nan;
%			dz(i,j)=nan;
%			dz_correct(i,j)=nan;
%			correctV(i,j)=nan;
%			gmtdz(i,j)=nan;
%			gmtdz_correct(i,j)=nan;
%			gmtcorrectV(i,j)=nan;
			%ampmap(i,j)=nan;
			GV(j,i)=NaN;
		end
	end
end


%% Calculate the error of surface fitting
%[m n]=size(xi);
%
%	%find out the grids where the stations are
%for k=1:length(coor)
%	staj=find(abs(xi(1,:)-coor(k,1))<gridsize/2);
%	if length(staj)>1
%		staj=staj(1);
%	end
%	stai=find(abs(yi(:,1)-coor(k,2))<gridsize/2);
%	if length(stai)>1
%		stai=stai(1);
%	end
%	stamatrix(stai,staj)=1;
%	staerr(k)=z(stai,staj)-coor(k,3);
%	gmtstaerr(k)=gmtz(stai,staj)-coor(k,3);
%end
%
%	% check the azi distribution of nearby station
%	% azicheck array is difined by:
%	% 2 | 1
%	% --+--
%	% 3 | 4
%	%
%stad=2; % check range by degree
%[m n]=size(xi);
%relmatrix=zeros(m,n);
%for i=1:m
%	for j=1:n
%		azicheck(1:4)=0;
%		if ~isnan(dz(i,j))
%			staj=find(abs(coor(:,1)-xi(i,j))<stad);
%			stai=find(abs(coor(staj,2)-yi(i,j))<stad/cosd(xi(i,j)));
%			staj=staj(stai);
%			stai=find(dist(coor(staj,1),coor(staj,2),xi(i,j),yi(i,j))<stad*111);
%			stai=staj(stai);
%			if ~isempty(stai)
%				for k=1:length(stai)
%					dla=coor(stai(k),1)-xi(i,j);
%					dlo=coor(stai(k),2)-yi(i,j);
%					dlo=dlo*cosd(xi(i,j));
%					if dla>0 && dla<stad
%						if dlo>0 && dlo<stad
%							azicheck(1)=1;
%						elseif dlo>-stad
%							azicheck(2)=1;
%						end
%					elseif dla>-stad
%						if dlo>0 && dlo<stad
%							azicheck(4)=1;
%						elseif dlo>-stad
%							azicheck(3)=1;
%						end
%					end
%				end
%%				if sum(azicheck(1:4))>3.9
%%					relmatrix(i,j)=1;
%%				end
%				relmatrix(i,j)=sum(azicheck(1:4));
%			end
%		end
%	end
%end
%
%
%for i=1:m
%	for j=1:n
%		if relmatrix(i,j)<3.5 && ~isnan(dz(i,j))
%			aspectch(i,j)=nan;
%			dz(i,j)=nan;
%			dz_correct(i,j)=nan;
%			correctV(i,j)=nan;
%			%ampmap(i,j)=nan;
%		end
%	end
%end
%
% save the workspace
filename = sprintf('%s_%1d.mat',event,period);
save(filename)


%if isonefigure
%	drawevent(event,period);
%	figure(2)
%	clf
%		ax = usamap(lalim, lolim);
%		set(ax, 'Visible', 'off')
%		states = shaperead('usastatehi', 'UseGeoCoords', true);
%		geoshow(ax, states, 'FaceColor', [0.5 0.5 1])
%		surfacem(xi,yi,gmtdz_correct);
%		contourm(xi,yi,gmtz,30,'k');
%		geoshow(ax, states, 'FaceColor', 'none')
%		seiscolormap
%		load seiscmap
%		colormap(seiscmap);
%		%quiverm(bigxi,bigyi,u,v);
%		caxis(phvrange(period+1,:));
%		colorbar;
%		plotm(coor(:,1),coor(:,2),'rv');
%else
%
%	figure(1);
%	clf
%		ax = usamap(lalim, lolim);
%		set(ax, 'Visible', 'off')
%		states = shaperead('usastatehi', 'UseGeoCoords', true);
%		geoshow(ax, states, 'FaceColor', [0.5 0.5 1])
%		surfacem(xi,yi,dz);
%		contourm(xi,yi,z,30,'k');
%		geoshow(ax, states, 'FaceColor', 'none')
%		seiscolormap
%		load seiscmap
%		colormap(seiscmap);
%		%quiverm(bigxi,bigyi,u,v);
%		caxis(phvrange(period+1,:));
%		colorbar;
%		plotm(coor(:,1),coor(:,2),'rv');
%		filename = sprintf('%s_%1d_grdmap',event,period);
%		title(filename,'Interpreter','none');
%		print('-depsc',filename);
%
%	% Calculate the aspect change compare to the event back azimuth
%	figure(2);
%	clf
%		ax = usamap(lalim, lolim);
%		set(ax, 'Visible', 'off')
%		states = shaperead('usastatehi', 'UseGeoCoords', true);
%		geoshow(ax, states, 'FaceColor', [0.5 0.5 1])
%		surfacem(xi,yi,aspectch);
%		%contourm(xi,yi,z,30,'k');
%		geoshow(ax, states, 'FaceColor', 'none')
%		%caxis([2 5]);
%		colorbar;
%		plotm(coor(:,1),coor(:,2),'rv');
%		filename=sprintf('Surface Surface Gradient Direction event:%s PID:%d',event,period);
%		title(filename,'Interpreter','none');
%		filename = sprintf('%s_%1d_azimap',event,period);
%		print('-depsc',filename);
%
%	figure(3);
%	clf
%		ax = usamap(lalim, lolim);
%		set(ax, 'Visible', 'off')
%		states = shaperead('usastatehi', 'UseGeoCoords', true);
%		geoshow(ax, states, 'FaceColor', [0.5 0.5 1])
%		surfacem(xi,yi,ampmap);
%		%caxis([2 5]);
%		colorbar;
%		plotm(coor(:,1),coor(:,2),'rv');
%		filename=sprintf('Amplitude Map event:%s PID:%d',event,period);
%		title(filename,'Interpreter','none');
%		filename = sprintf('%s_%1d_ampmap',event,period);
%		print('-depsc',filename);
%
%	figure(4)
%	clf
%		ax = usamap(lalim, lolim);
%		set(ax, 'Visible', 'off')
%		states = shaperead('usastatehi', 'UseGeoCoords', true);
%		surfacem(xi,yi,correctV);
%		%contourm(xi,yi,z,30,'k');
%		geoshow(ax, states, 'FaceColor', 'none')
%		sc=nanmean(nanstd(correctV));
%		caxis([-4*sc 4*sc]);
%		seiscolormap
%		load seiscmap
%		colormap(seiscmap);
%		colorbar;
%		%plotm(coor(:,1),coor(:,2),'rv');
%		filename=sprintf('Amplitude Correction Map event:%s PID:%d',event,period);
%		title(filename,'Interpreter','none');
%		filename = sprintf('%s_%1d_ampcorrection',event,period);
%		print('-depsc',filename);
%
%	figure(5);
%	clf
%		ax = usamap(lalim, lolim);
%		set(ax, 'Visible', 'off')
%		states = shaperead('usastatehi', 'UseGeoCoords', true);
%		geoshow(ax, states, 'FaceColor', [0.5 0.5 1])
%		surfacem(xi,yi,dz_correct);
%		contourm(xi,yi,z,30,'k');
%		geoshow(ax, states, 'FaceColor', 'none')
%		seiscolormap
%		load seiscmap
%		colormap(seiscmap);
%		%quiverm(bigxi,bigyi,u,v);
%		caxis(phvrange(period+1,:));
%		colorbar;
%		plotm(coor(:,1),coor(:,2),'rv');
%		filename = sprintf('%s_%1d_grdmap_correct',event,period);
%		title(filename,'Interpreter','none');
%		print('-depsc',filename);
%
%end % End of figure plot if


end	% end of period loop
event=fgetl(fp);

end % end of event loop
