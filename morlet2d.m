function mo=morlet2d(a,theta)
% This function is used to generate a 2D morlet wavelet matrix.
% Ge Jin, jinwar@gmail.com

	N=300;
	k_0=5.336;
% x-y domain for physical space wavelet
	x = [-N/2:N/2 - 1]*(1/N); % periodic domain size = 1, from -1/2 -> 1/2
	xs = x/a; % dilation of wavelet
	[x,y] = meshgrid(xs); % 2D grid of x’s and y’s (x,y)
% rotate the coordinates
	xd=(x*cos(theta)+y*sin(theta));
	yd=(y*cos(theta)-x*sin(theta));
	r2=xd.^2+yd.^2;
% generate a 2D complex Morlet wavelet in physical space
	% corrected wavelet
	mo = ( exp(sqrt(-1)*k_0*xd) - exp(-k_0^2/2) ) .* exp(-r2/2);
	mo=(1/a)*mo; % renormalize due to dilation in **2D**
	
