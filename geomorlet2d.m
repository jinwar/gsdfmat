function mo=geomorlet2d(xi,yi,lat0,lon0,lamda,theta)
% This function is used to generate a 2D morlet wavelet matrix for a geo-coor.
% xi yi are coors, lat0, lon0 are the center of the filter, lamda is wavelength, theta is azimuth of plain wave
% Ge Jin, jinwar@gmail.com

	l=3*lamda;
	a=0.1329*l;
	k_0=2*pi./lamda.*a;
	theta=deg2rad(theta);
	%theta=deg2rad(theta-90);

	[m,n]=size(xi);

	x=xi-lat0;
	y=yi-lon0;

	x=x.*111.1949;
	for i=1:m
		for j=1:n
			y(i,j)=y(i,j)*cos(deg2rad(xi(i,j)))*111.1949;
		end
	end

	r2=dist(xi,yi,lat0,lon0);

% x-y domain for physical space wavelet
%	x = [-N/2:N/2 - 1]*(1/N); % periodic domain size = 1, from -1/2 -> 1/2
%	xs = x/a; % dilation of wavelet
%	[x,y] = meshgrid(xs); % 2D grid of x’s and y’s (x,y)

	x=x/a;y=y/a;
	r2=r2/a;

% rotate the coordinates
	xd=(x*cos(theta)+y*sin(theta));
	yd=(y*cos(theta)-x*sin(theta));

% generate a 2D complex Morlet wavelet in physical space
	% corrected wavelet
	mo = ( exp(sqrt(-1)*k_0*xd) - exp(-k_0^2/2) ) .* exp(-r2/2);
	mo=(1/a)*mo; % renormalize due to dilation in **2D**
	
