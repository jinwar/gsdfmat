function cwt=CWT_geomorlet2d(xi,yi,z,lamda,theta)
% Function to make the 2d morlet transform for a geographic map
% calculate in physical domain, extremely slow!
% Ge Jin, jinwar@gmail.com
%

[m,n]=size(z);
cwt=zeros(m,n);
areamat=zeros(m,n);
z=z-nanmean(nanmean(z));

% Build up area matrix
for i=1:m-1
	for j=1:n-1
		dlat=abs(xi(i+1,j+1)-xi(i,j));
		dlon=abs(yi(i+1,j+1)-yi(i,j));
		areamat(i,j)=dlat*dlon*cos(xi(i,j))*111.19^2;
	end
end
areamat(m,:)=areamat(m-1,:);
areamat(:,n)=areamat(:,n-1);

sumarea=sum(sum(areamat));

for i=1:m
	i
	for j=1:n
		if isnan(z(i,j))
			cwt(i,j)=NaN;
		else
			mo=geomorlet2d(xi,yi,xi(i,j),yi(i,j),lamda,theta);
			mo=mo.*z;
			mo=mo.*areamat;
			cwt(i,j)=nansum(nansum(mo));
		end
	end
end

cwt=cwt./sumarea;

