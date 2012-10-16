function checksta(eventname,period)

checkcs=1;
checkwaveform=0;

filename=sprintf('%s_%1d.mat',eventname,period);
load(filename);

figure(98);
clf
	ax = usamap(lalim, lolim);
	set(ax, 'Visible', 'off')
	states = shaperead('usastatehi', 'UseGeoCoords', true);
	geoshow(ax, states, 'FaceColor', [0.5 0.5 1])
	surfacem(xi,yi,dz_correct);
	contourm(xi,yi,z,20,'k');
	geoshow(ax, states, 'FaceColor', 'none')
	seiscolormap
	load seiscmap
	colormap(seiscmap);
	%quiverm(bigxi,bigyi,u,v);
% 	caxis(phvrange(period+1,:));
	colorbar;
	plotm(coor(:,1),coor(:,2),'bv','markersize',5);
	filename = sprintf('%s_%1d_grdmap_correct',event,period);
	title(filename,'Interpreter','none');


	[mla mlo]=inputm(1);

	ddist=distance(stadata(:,2),stadata(:,3),mla,mlo);

n=find(ddist==min(ddist));
staindex=stadata(n,1);
stadata(n,2:3)
stla(1)=stadata(n,2);
stlo(1)=stadata(n,3)+360;
staepidist=distance(stadata(n,2),stadata(n,3),evla,evlo);

stemp=sprintf('cat %s.sta | awk ''{if ($2==%d) print $1 }'' > tempfile ',eventname,staindex);
system(stemp);
staname=textread('tempfile','%s')

plotm(stadata(n,2),stadata(n,3)+360,'rv');
for i=1:length(newcsdata)
    if newcsdata(i,1)==staindex 
        err=newcsdata(i,3)-(tnet(floor(newcsdata(i,1))+1)-tnet(floor(newcsdata(i,2))+1));
        stla(2)=stadata(newcsdata(i,2)+1,2);
        stlo(2)=stadata(newcsdata(i,2)+1,3)+360;
        if err<-0.2
            plotm(stla,stlo,'r');
        elseif err>0.2
            plotm(stla,stlo,'b');
        else
            plotm(stla,stlo,'g');
        end
    elseif newcsdata(i,2)==staindex 
        err=newcsdata(i,3)-(tnet(floor(newcsdata(i,1))+1)-tnet(floor(newcsdata(i,2))+1));
        stla(2)=stadata(newcsdata(i,1)+1,2);
        stlo(2)=stadata(newcsdata(i,1)+1,3)+360;
        if err<-0.2
            plotm(stla,stlo,'r');
        elseif err>0.2
            plotm(stla,stlo,'b');
        else
            plotm(stla,stlo,'g');
        end
    end
end	


staids(1,1)=staindex;
staids(1,2)=staepidist;
staids(1,3)=stla(1);
staids(1,4)=stlo(1);
staids(1,5)=0;
staids(1,6)=stadata(staindex+1,period+4);
stan=1;
if checkcs
	figure(99)
    clf
	hold on
	for i=1:length(csdata)
		if csdata(i,1)==staindex 
			err=csdata(i,3)-(tnet(floor(csdata(i,1))+1)-tnet(floor(csdata(i,2))+1));
            stla=stadata(csdata(i,2)+1,2);
            stlo=stadata(csdata(i,2)+1,3);
            epidist=distance(stla,stlo,evla,evlo);
			plot(epidist-staepidist,err,'x');
		elseif csdata(i,2)==staindex 
			err=csdata(i,3)-(tnet(floor(csdata(i,1))+1)-tnet(floor(csdata(i,2))+1));
            stla=stadata(csdata(i,1)+1,2);
            stlo=stadata(csdata(i,1)+1,3);
            epidist=distance(stla,stlo,evla,evlo);
			plot(epidist-staepidist,err,'x');
		end
    end	
    for i=1:length(newcsdata)
		if newcsdata(i,1)==staindex 
			err=newcsdata(i,3)-(tnet(floor(newcsdata(i,1))+1)-tnet(floor(newcsdata(i,2))+1));
            stla=stadata(newcsdata(i,2)+1,2);
            stlo=stadata(newcsdata(i,2)+1,3);
            epidist=distance(stla,stlo,evla,evlo);
            plot(epidist-staepidist,err,'rx');
			stan=stan+1;
			staids(stan,1)=newcsdata(i,2);
			staids(stan,2)=epidist;
			staids(stan,3)=stla;
			staids(stan,4)=stlo;
			staids(stan,5)=err;
			staids(stan,6)=stadata(newcsdata(i,2)+1,period+4);

		elseif newcsdata(i,2)==staindex 
			err=newcsdata(i,3)-(tnet(floor(newcsdata(i,1))+1)-tnet(floor(newcsdata(i,2))+1));
            stla=stadata(newcsdata(i,1)+1,2);
            stlo=stadata(newcsdata(i,1)+1,3);
            epidist=distance(stla,stlo,evla,evlo);
            plot(epidist-staepidist,err,'rx');
			stan=stan+1;
			staids(stan,1)=newcsdata(i,1);
			staids(stan,2)=epidist;
			staids(stan,3)=stla;
			staids(stan,4)=stlo;
			staids(stan,5)=err;
			staids(stan,6)=stadata(newcsdata(i,1)+1,period+4);
		end
	end	
end
figure(97)
clf
plot(staids(:,2),staids(:,6),'x')

% write the sac scripts to read in the waveforms and plot them

staids=sortrows(staids,2);

fp=fopen('sacmacro.csh','w');

fprintf(fp,'sac<<!\n');
stemp=sprintf('cat %s.sta | awk ''{if ($2==%d) print $1 }'' > tempfile ',eventname,staids(1,1));
system(stemp);
staname=textread('tempfile','%s');
fprintf(fp,'r more %s/*%s*LHZ.sac.w\n',eventname,char(staname));
for i=1:length(staids)
	stemp=sprintf('cat %s.sta | awk ''{if ($2==%d) print $1 }'' > tempfile ',eventname,staids(i));
    system(stemp);
    staname=textread('tempfile','%s');
	fprintf(fp,'r more %s/*%s*LHZ.sac_%1d.bp\n',eventname,char(staname),period);
end
fprintf(fp,'r more %s/*%s*LHZ.sac.w\n',eventname,char(staname));
fprintf(fp,'ppk\n');
fprintf(fp,'q\n');
fprintf(fp,'!\n');

fclose(fp);

if checkwaveform
    system('csh sacmacro.csh');
end


end
