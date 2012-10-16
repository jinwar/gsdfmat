
	G=phaseg;
	[xi yi]=ndgrid(xnode,ynode);

	Nx=length(xnode);
	Ny=length(ynode);
%	nray=size(G,1);
	GVx=zeros(Nx,Ny);
	GVy=zeros(Nx,Ny);

	for i=1:Nx
		for j=1:Ny
			n=Ny*(i-1)+j;
			GVx(i,j)= G(2*n-1);
			GVy(i,j)= G(2*n);
		end
	end

	GV=(GVx.^2+GVy.^2).^-.5;
    lalim=[min(xnode) max(xnode)];
    lolim=[min(ynode) max(ynode)];
	

	err = mat*phaseg - dt;
	err = W*err;

	inderr=find(abs(err)>2*std(err));

	figure(4)
	clf
	hold on
		ax = usamap(lalim, lolim);
		set(ax, 'Visible', 'off')
		states = shaperead('usastatehi', 'UseGeoCoords', true);
		surfacem(xi,yi,GV)
%		plotm(stadata(:,2),stadata(:,3),'rv');
		%contourm(xi,yi,z,30,'k');
		geoshow(ax, states, 'FaceColor', 'none')
		seiscolormap
		load seiscmap
		colormap(seiscmap);
		caxis([3.5 4.15])

		colorbar

	for i=1:length(inderr)
		sta1id=csdata(inderr(i),1)+1;
		sta2id=csdata(inderr(i),2)+1;
		errla(1)=stadata(sta1id,2);
		errla(2)=stadata(sta2id,2);
		errlo(1)=stadata(sta1id,3);
		errlo(2)=stadata(sta2id,3);
		plotm(errla,errlo,'k')
	end

	figure(5)
	hist(err,100);
