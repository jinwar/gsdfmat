function plotkernel(G,rays, xnode,ynode)

	[xi yi]=ndgrid(xnode,ynode);

	Nx=length(xnode);
	Ny=length(ynode);
	nray=size(G,1);
	GVx=zeros(Nx,Ny);
	GVy=zeros(Nx,Ny);
	karray=1:30:nray;

	sumG=sum(abs(G),1);

    for i=1:Nx
        for j=1:Ny
            n=Ny*(i-1)+j;
            GVx(i,j)=sumG(2*n-1);
            GVy(i,j)=sumG(2*n);
        end
    end

	GV=(GVx.^2+GVy.^2).^.5;
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
		%contourm(xi,yi,z,30,'k');
		geoshow(ax, states, 'FaceColor', 'none')
		colorbar

	figure(3)
	clf
		ax = usamap(lalim, lolim);
		set(ax, 'Visible', 'off')
		states = shaperead('usastatehi', 'UseGeoCoords', true);
		surfacem(xi,yi,GV)
		%contourm(xi,yi,z,30,'k');
		geoshow(ax, states, 'FaceColor', 'none')
		for k=karray
			lats(1)=rays(k,1);
			lats(2)=rays(k,3);
			longs(1)=rays(k,2);
			longs(2)=rays(k,4);
			plotm(lats,longs,'r');
			plotm(lats,longs,'rv');
        end
		colorbar
