function cwt=FFTCWT_geomorlet2d(xi,yi,z,lamda,theta)
% Function to make the 2d morlet transform for a geographic map
% Ge Jin, jinwar@gmail.com
%

[m,n]=size(z);
cwt=zeros(m,n);
areamat=zeros(m,n);
z=z-nanmean(nanmean(z));

for i=1:m
	for j=1:n
		if isnan(z(i,j))
			z(i,j)=0;
		end
	end
end

lat0=mean(mean(xi));
lon0=mean(mean(yi));

mo=geomorlet2d(xi,yi,lat0,lon0,lamda,theta);

fftmo=fft2(mo);
fftz=fft2(z);


%fftmo=conj(fftmo).*fftz;
fftmo=conj(fftmo).*(fftz);

%fftmo=fftmo.*exp(-sqrt(-1)*pi/4);

cwt=fftshift(ifft2(fftmo));

