function output=stainform(pid)

	IsFigure=0

	IsMap=1

	lalim=[28 50.5];
	lolim=[-115 -95];
	lolim=lolim+360;
	xnode=lalim(1):0.1:lalim(2);
	ynode=lolim(1):0.1:lolim(2);


	system('ls stainform/*_2.phv | awk -F_ ''{print $1}'' > stationlist');

	fp=fopen('stationlist','r');

	staname=fgetl(fp);

	stan=0;

	while ischar(staname)
		filename=sprintf('%s_%1d.phv',staname,pid)
		data=load(filename);
		[m n]=size(data);
		evn=0;
		clear evdata;
		for i=1:m
			if data(i,10)>0
				evn=evn+1;
				% evdata(:,1)=azi; evdata(:,2)=v; evdata(:,3)=std(v);
				evdata(evn,1)=azimuth(data(i,3),data(i,4),data(i,1),data(i,2));
				evdata(evn,2)=data(i,5);
				if data(i,6)>0
					evdata(evn,3)=data(i,6);
				else
					evdata(evn,3)=0.1;
                end
			end
		end
		if evn >0
			stan=stan+1;
			stainf(stan,1)=data(1,3);
			stainf(stan,2)=data(1,4);
			stainf(stan,3)=mean(evdata(:,2));
			%stainf(stan,3)=weightavg(evdata);
			para=[stainf(stan,3) 0 0]; % No anisotropy
			%para=gridsearch(evdata);
			
			if stainf(stan,3) < 3.7
				IsFigure =0;
			end
			if IsFigure
				x=0:359;
				v=anisoV(para,x);
				figure(1)
				clf;
				hold on;
				errorbar(evdata(:,1),evdata(:,2),evdata(:,3),'x');
				plot(x,v,'r');
				keyboard
			end
		end
		staname=fgetl(fp);
	end % end of station

	if IsMap

		coor=stainf(:,1:2);
		coor(:,2)=coor(:,2)+360;
		[z,xi,yi]=gridfit(coor(:,1),coor(:,2),stainf(:,3),xnode,ynode,...
							'smooth',2,'solver','normal');
		interpz=griddata(coor(:,1),coor(:,2),stainf(:,3),xi,yi,'cubic');
		[m,n]=size(interpz);
		for i=1:m
			for j=1:n
				if isnan(interpz(i,j))
					z(i,j)=nan;
				end
			end
		end

		figure(1);
		clf
			usamap(lalim,lolim);
			[topoZ, refvec] = etopo('ETOPO5.DAT',1,lalim,lolim);
			geoshow(topoZ, refvec,'DisplayType','texturemap');
			colormap(demcmap(topoZ));
			plotm(coor(:,1),coor(:,2),'rv');
			title('Topography');

		figure(2);
		clf
			ax = usamap(lalim, lolim);
			set(ax, 'Visible', 'off')
			states = shaperead('usastatehi', 'UseGeoCoords', true);
			geoshow(ax, states, 'FaceColor', [0.5 0.5 1])
			surfacem(xi,yi,z);
			seiscolormap
			load seiscmap
			colormap(seiscmap);
			colorbar;
			plotm(coor(:,1),coor(:,2),'rv');
			filename=sprintf('Station Phase Velocity (gridfit), period:%d',pid)
			title(filename);
			filename=sprintf('staphvmap_%1d.png',pid)
			print('-dpng',filename);

	end

end % end of function

function para=gridsearch(evdata)
	meanv=weightavg(evdata);
	v0=0.9:0.02:1.1;
	dv=0:0.02:0.1;
	theta=0:5:180;
	v0=v0*meanv;
	dv=dv*meanv;

	minerr=inf;
	minindex=[0 0 0];
	for i=1:length(v0)
		for j=1:length(dv)
			for k=1:length(theta)
				err(i,j,k) = errfun([v0(i) dv(j) theta(k)],evdata);
				if minerr > err(i,j,k)
					minerr = err(i,j,k);
					minindex=[i j k];
				end
			end
		end
	end
	i=minindex(1); j=minindex(2); k=minindex(3);
	para(1)=v0(i); para(2)=dv(j); para(3)=theta(k);
	if 0
		figure(2)
		clf
		[x,y]=meshgrid(theta,dv);
		for j=1:length(dv)
			for k=1:length(theta)
				z(j,k)=err(i,j,k);
			end
		end
		contour(x,y,z);
		colorbar;
	end
	para
end

function v0=weightavg(evdata)
	[m n]=size(evdata);
	sumvar=0;
	sumv0=0;
	for i=1:m
		sumv0=sumv0+evdata(i,2)/(evdata(i,3)^2);
		sumvar=sumvar+1/(evdata(i,3)^2);
	end
	v0=sumv0/sumvar;
end
% The anisotropic velocity can be presented by 
% v(azi)=V0+dv*cos(azi - theta)
% v0=para(1); dv=para(2); theta=para(3);
function err=errfun(para,evdata)
	v0=para(1); dv=para(2); theta=para(3);
	[m,n]=size(evdata);
	sumerr=0;
	sumvar=0;
	for i=1:m
		e=evdata(i,2)-(v0+dv*cosd(2*evdata(i,1)-theta));
		sumerr=sumerr+(e/evdata(i,3)).^2;
		sumvar=sumvar+(1/evdata(i,3)).^2;
	end
	err=sumerr/sumvar;
end

function v=anisoV(para,azi)
	v0=para(1); dv=para(2); theta=para(3);
	v=v0+dv*cosd(2*azi-theta);
end
