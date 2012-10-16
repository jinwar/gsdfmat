event='201010200409'

system('ls 201010200409/*2.stadt > stadtlist');

fp=fopen('stadtlist','r')

stadtfile=fgetl(fp);
netnum=1;

sta1.la=43.999;
sta1.lo=-112.485;
sta2.la=42.761;
sta2.lo=-100.318;

while ischar(stadtfile)
	stadtdata=load(stadtfile);
	[m n]=size(stadtdata);
	for i=1:m
		if distance(stadtdata(i,1),stadtdata(i,2),sta1.la,sta1.lo) < 0.05 
			sta1.dt=stadtdata(i,3);
		end
		if distance(stadtdata(i,1),stadtdata(i,2),sta2.la,sta2.lo) < 0.05 
			sta2.dt=stadtdata(i,3);
		end
	end
	dt(netnum)=sta1.dt-sta2.dt;
	disp(dt(netnum))
	netnum=netnum+1;
	stadtfile=fgetl(fp)
end
