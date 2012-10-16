% This program is used to exclude close event that might generate interference
% with each other

clear
evfp=fopen('alleventlist','r');

eventid=fgetl(evfp);
eventnum=0;
groupv=3.3;
stla=40;
stlo=-110;
outfile='eventlist';

while ischar(eventid)
	eventnum=eventnum+1;
	[y m d h min]=strread(eventid,'%4d%2d%2d%2d%2d');
	eventinfo(eventnum).id=eventid;
	eventinfo(eventnum).otime=datenum(y,m,d,h,min,0)*24*3600;
	stemp=['ls ',eventid,'/*LHZ*sac | head -n1 > sactemplist'];
	system(stemp);
	sacfp=fopen('sactemplist','r');
	sacfilename=fgetl(sacfp);
	fclose(sacfp);
	sac=readsac(sacfilename);
	eventinfo(eventnum).evla=sac.EVLA;
	eventinfo(eventnum).evlo=sac.EVLO;
	dist=deg2km(distance(sac.EVLA,sac.EVLO,stla,stlo));
	eventinfo(eventnum).atime=eventinfo(eventnum).otime + dist./groupv;
	eventinfo(eventnum).isgood=1;
	eventid=fgetl(evfp)
end

fclose(evfp);

for i=2:eventnum-1;
	if eventinfo(i-1).atime > eventinfo(i).atime - 2000
		eventinfo(i-1).isgood=0;
		eventinfo(i).isgood=0;
	end
	if eventinfo(i+1).atime < eventinfo(i).atime + 2000
		eventinfo(i+1).isgood=0;
		eventinfo(i).isgood=0;
	end
end

evfp=fopen(outfile,'w')
for i=1:eventnum
	if eventinfo(i).isgood==1
	fprintf(evfp,'%s\n',eventinfoi(i).id);
	end
end
fclose(evfp);

	
