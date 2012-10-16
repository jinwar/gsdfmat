function smap=smoothmap_old(xi,yi,map,D)
% This function is used to make a gaussian smoothing on a grd map

% D is Gaussian average distance

[m n]=size(xi);

midi=floor(m/2);
midj=floor(n/2);

dx=vdist(xi(midi,midj),yi(midi,midj),xi(midi+1,midj),yi(midi+1,midj))/1e3;
dy=vdist(xi(midi,midj),yi(midi,midj),xi(midi,midj+1),yi(midi,midj+1))/1e3;

Nx=floor(D/dx);
Ny=floor(D/dy);
w=zeros(2*Nx+1,2*Ny+1);
for i=-Nx:Nx
	for j=-Ny:Ny
		if i==0 & j==0
			w(i+Nx+1,j+Ny+1)=1;
		else
			dist=vdist(xi(midi,midj),yi(midi,midj),xi(midi+i,midj+j),yi(midi+i,midj+j))/1e3;
			w(i+Nx+1,j+Ny+1)=exp(-dist.^2./D.^2/2);
		end
	end
end

sumall=sum(sum(w));
w=w./sumall;

smap=conv2(map,w,'same');
