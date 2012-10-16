function plotphaseg(G,stadata,xnode,ynode)

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
	
	figure(1)
	clf
		ax = usamap(lalim, lolim);
		set(ax, 'Visible', 'off')
		states = shaperead('usastatehi', 'UseGeoCoords', true);
		surfacem(xi,yi,GVx)
		%contourm(xi,yi,z,30,'k');
		geoshow(ax, states, 'FaceColor', 'none')
		colorbar

	figure(2)
	clf
		ax = usamap(lalim, lolim);
		set(ax, 'Visible', 'off')
		states = shaperead('usastatehi', 'UseGeoCoords', true);
		surfacem(xi,yi,GVy)
		plotm(stadata(:,2),stadata(:,3),'rv');
		%contourm(xi,yi,z,30,'k');
		geoshow(ax, states, 'FaceColor', 'none')
		colorbar

	figure(3)
	clf
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
