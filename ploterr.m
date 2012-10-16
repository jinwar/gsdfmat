
goodnum=0;
for i=1:length(out)
	%if (out(i,12)<4 && out(i,14)>0.5 && out(i,6)==150 )
	%if (out(i,12)<4 && out(i,14)>0.5 )
	%if (out(i,12)<200 && out(i,13)>2)
	if 1
		goodnum=goodnum+1;
		gooddata(goodnum,1)=out(i,7);
		gooddata(goodnum,2)=out(i,5);
		gooddata(goodnum,3)=out(i,6);
		%gooddata(goodnum,2)=distance(out(i,1),out(i,2),out(i,3),out(i,4))*111;
%		if abs(out(i,5))>1
%			gooddata(goodnum,2)=out(i,7);
%		else
%			gooddata(goodnum,2)=out(i,7);
%		end
	end
end
goodnum/length(out)

coor=out(:,1:4);

lalim=[min(coor(:,1))-1 max(coor(:,1))+1];
lolim=[min(coor(:,2))-1 max(coor(:,2))+1];
lolim=lolim+360;

if isplotmap
	figure(1);
	clf
	usamap(lalim,lolim);
	[Z, refvec] = etopo('ETOPO5.DAT',5,lalim,lolim);
	geoshow(Z, refvec,'DisplayType','texturemap');
	colormap(demcmap(Z));

	hold on
	for i=1:length(coor)
		x(1)=coor(i,1);
		x(2)=coor(i,3);
		y(1)=coor(i,2);
		y(2)=coor(i,4);
		linem(x,y);
	end
	plotm(sta(:,1),sta(:,2),'rv');
end

periods=[20 30 40 60 100 150];
%plot(out(:,10),out(:,6),'x');
for i=1:length(periods)
	figure(i)
	clf
	hold on
	for j=1:length(gooddata)
		if gooddata(j,3)==periods(i)
			plot(gooddata(j,1),gooddata(j,2),'x');
		end
	end
end
for i=1:length(periods)
	figure(i)
	title(sprintf('period: %d s', periods(i)))
	xlim([-20 20])
	ylim([-250 250])
	filename=sprintf('err_%d',i);
	print('-djpeg',filename);
end

%xlim([90 100]);
%ylim([-10 10])
%xlim([0 120]);
%print('-dpng','cferr')
