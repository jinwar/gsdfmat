% Matlab program to measure the phase and amplitude of surface wave using FTAN method descripted by
% Levshin, 1992
%
% Written by Ge Jin, jinwar@gmail.com
% LDEO, columbia University
% 2011-12
%
clear

eventslist='events'

minp=25;
maxp=100;
np=20;
for i=1:np
	periods(i)=minp+exp((log(maxp-minp+1))/(np-1)*(i-1))-1;
end
%periods=[25 32 40 50 66 83 100];
freqs=1./periods;
minalpha=1e3;
maxalpha=1e4;
for i=1:np
	alpha(i)=minalpha+exp((log(maxalpha-minalpha+1))/(np-1)*(i-1))-1;
end

fpev=fopen(eventslist,'r');
eventline=fgetl(fpev);
eventnum=0;

Isfigure=1;

while ischar(eventline) 
	[eventid v1 t1 v2 t2]=strread(eventline,'%s %f %f %f %f\n');

	% Get sac file list
	stemp=sprintf('ls %s/*.TA.*LHZ.sac > tempstalist',char(eventid));
	system(stemp);
	fpsta=fopen('tempstalist','r');
	stafile=fgetl(fpsta);
	stanum=1;
	clear stainfo stadata stadist;

	while ischar(stafile)
		sac=readsac(stafile);
		epidist=sac.DIST;
		winbegin=sac.DIST/v1+t1;
		winend=sac.DIST/v2+t2;

		sacdata=sac.DATA1';
		% set time axis
		sact=sac.B:sac.DELTA:sac.B+(sac.NPTS-1)*sac.DELTA;

		% window the data
		sacwindata=sacdata;
		for i=1:length(sact)
			if sact(i) < winbegin || sact(i) > winend
				sacwindata(i)=0;
			end
		end

		%fft and frequency axis
		sacfft=fft(sacwindata);
		faxis=0:1/(sact(end)-sact(1)):1/sac.DELTA/2;
		sacffthalf=sacfft(1:length(faxis));

		% Generate narrow band gaussian filtered data
		clear gausf nband
		for i=1:length(freqs)
			gausf(i,:)=exp(-alpha(i)*(faxis-freqs(i)).^2./2./freqs(i));
			nbpfft=sacffthalf.*gausf(i,:);
			nbpfft(end+1:length(sacfft))=0;
			nband(i,:)=ifft(nbpfft);
			maxamp(i)=find(abs(nband(i,:))==max(abs(nband(i,:))));
		end

		if Isfigure
			figure(1)
			clf;
			[yi xi]=meshgrid(sac.DIST./sact,periods);
			contour(xi,yi,abs(nband));
			hold on
			for i=1:length(freqs)
				maxampx=periods(i);
				maxampy=sac.DIST./sact(maxamp(i));
				plot(maxampx,maxampy,'x');
			end
			hold off
			set(gca,'xscale','log');
			set(gca,'xtick',[25 30 40 60 100]);
			shading flat
			ylim([3 4.5])
			xlim([minp maxp])
			pause
		end

		stafile=fgetl(fpsta);
		stanum=stanum+1;
	end % end of loop stations

	fclose(fpsta);
	eventline=fgetl(fpev);
end % end of event loop
