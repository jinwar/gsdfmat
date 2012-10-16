clear;

data=load('logsumfile');

evla=data(:,1);
evlo=data(:,2);
stla=36;
stlo=-120;

[dist,az]=distance(stla,stlo,evla,evlo);

periods=[20 30 40 60 100 150];
for j=1:length(periods)
	figure(j)
	clf
	hold on
	for i=1:length(data)
		if data(i,2+j*2)>100
			plot(az(i),data(i,1+2*j),'.');
		end
	end
	xlim([0 360])
	stemp=sprintf('Period = %d',periods(j));
	title(stemp)
	stemp=sprintf('phvaz_%1d',j);
	print('-dpng',stemp);
end

