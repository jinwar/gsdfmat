clear

for k=1:200
n=10;
v=4;
noise=5;

for i=1:n
	x(i)=(rand()-0.5)*500;
	y(i)=x(i)/v+ (rand()-0.5)*noise*2;
end
data(:,1)=x;
data(:,2)=y;

save data.txt data -ASCII

fit0(k,:)=polyfit(x,y,1);

system('./fitphasev data.txt');

fp=fopen('fitphasev.out','r');
fit1(k,:)=fscanf(fp,'%f %f %f\n');
fclose(fp);

end
v0=1./fit0(:,1);
v1=1./fit1(:,1);

disp('with constant')
mean(v0)
std(v0)

disp('without constant')
mean(v1)
std(v1)

