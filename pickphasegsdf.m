% This program is used to help to decide the window to generate isolate filter in the gsdf code.
clear

% parameter need to be set
eventslist='eventlist';
lalim=[30 50];
lolim=[-110 -90];
plotstanum=6;
periods=[25 40 60 100];
prefilter=[15 150];
groupv=[5 2];
refgv=[5 4 3 2];

fpev=fopen(eventslist,'r');
event=fgetl(fpev);
eventnum=0;
sleep=0;

while ischar(event) && sleep==0

	eventnum=eventnum+1;
	eventinfo(eventnum).name=event;
	stemp=sprintf('ls %s/*.TA.*LHZ.sac > tempstalist',event);
	system(stemp);
	fpsta=fopen('tempstalist','r');
	stafile=fgetl(fpsta);
	stanum=1;
	clear stainfo stadata stadist;

	while ischar(stafile)
		% Read in sac file
		sachdr=readsac(stafile);

		if sachdr.STLA > lalim(1) && sachdr.STLA < lalim(2)...   % test whether it's in the range
				&& sachdr.STLO > lolim(1) && sachdr.STLO < lolim(2)
			stainfo(stanum).dist=sachdr.DIST;
			stainfo(stanum).filename=stafile;
			stadist(stanum,1)=stanum;
			stadist(stanum,2)=sachdr.DIST;
			stanum=stanum+1;
		end
		stafile=fgetl(fpsta);
	end % Loop of stations

	fclose(fpsta)

	% sort the station by epicenter distance
	stadist_sort=sortrows(stadist,2);

	% select some stations to plot
	stanum=length(stadist);
	ddist=ceil((stadist_sort(end,2)-stadist_sort(1,2))/plotstanum);
	dist=stadist_sort(1,2):ddist:stadist_sort(end,2);  % dist is the array that stations being plotted.

	% Find and store the stations being plotted in the structure stadata
	for i=1:length(dist)
		[mind, stai]=min(abs(stadist(:,2)-dist(i)));
		stai=stadist(stai,1);
		sac=readsac(stainfo(stai).filename);
		stadata(i).filename=stainfo(stai).filename;
		stadata(i).data=sac.DATA1';
		stadata(i).data=stadata(i).data./max(abs((stadata(i).data))); % Normalize it
		stadata(i).dist=sac.DIST;
		stadata(i).stla=sac.STLA;
		stadata(i).stlo=sac.STLO;
		stadata(i).t=sac.B:sac.DELTA:sac.B+(sac.NPTS-1)*sac.DELTA;
	end

	% filter the waveforms
	fN=1/2/sac.DELTA;
	for i=1:length(dist)
		f1=1/prefilter(2);
		f2=1/prefilter(1);
		[b,a]=butter(2,[f1/fN, f2/fN]);
		stadata(i).prefdata=filtfilt(b,a,stadata(i).data);
		for j=1:length(periods)
			f1=1/periods(j)*.9;
			f2=1/periods(j)*1.1;
			[b,a]=butter(2,[f1/fN, f2/fN]);
			stadata(i).fdata(j,:)=filtfilt(b,a,stadata(i).data);
			stadata(i).fdata(j,:)=stadata(i).fdata(j,:)./max(abs(stadata(i).fdata(j,:)));
		end
	end
	
	% generate the envelope function
	for i=1:length(dist)
		for j=1:length(periods)
%			ind=find(diff(sign(diff(stadata(i).fdata(j,:))))==+2)+1;
%			stadata(i).envelop(j,:)=interp1(stadata(i).t(ind),stadata(i).fdata(j,ind),stadata(i).t,'linear');
			stadata(i).envelop(j,:)=abs(hilbert(stadata(i).fdata(j,:)));
		end
	end
	% Plot the original data
	figure(1)
	clf
		hold on
		amp=(dist(end)-dist(1))/plotstanum;
		for i=1:length(dist)
			plot(stadata(i).t, stadata(i).prefdata*amp+stadata(i).dist);
			for j=1:length(periods)
				plot(stadata(i).t, stadata(i).envelop(j,:)*amp/2+stadata(i).dist,'k');
			end
		end
		plot([dist(1)/groupv(1) dist(end)/groupv(1)],[dist(1) dist(end)],'r')
		plot([dist(1)/groupv(2) dist(end)/groupv(2)],[dist(1) dist(end)],'r')
		for k=1:length(refgv)
			plot([dist(1)/refgv(k) dist(end)/refgv(k)],[dist(1) dist(end)],'g--')
		end
		ylim([dist(1)-ddist dist(end)+ddist])
		xlim([dist(1)/groupv(1)-500 dist(end)/groupv(2)+1000])

	% Plot filtered data
	figure(2)
	clf
	for i=1:length(periods)
		subplot(1,length(periods),i)
			hold on
			amp=(dist(end)-dist(1))/plotstanum/2;
			for j=1:length(dist)
				plot(stadata(j).t, stadata(j).fdata(i,:)*amp+stadata(j).dist);
				plot(stadata(j).t, stadata(j).envelop(i,:)*amp+stadata(j).dist,'k');
			end
			plot([dist(1)/groupv(1) dist(end)/groupv(1)],[dist(1) dist(end)],'r')
			plot([dist(1)/groupv(2) dist(end)/groupv(2)],[dist(1) dist(end)],'r')
			for k=1:length(refgv)
				plot([dist(1)/refgv(k) dist(end)/refgv(k)],[dist(1) dist(end)],'g--')
			end
			ylim([dist(1)-ddist dist(end)+ddist])
			xlim([dist(1)/groupv(1)-500 dist(end)/groupv(2)+1000])
			title(sprintf('%d',periods(i)));
	end

	% Plot the world map and great circle path
	figure(3)
	clf
		evla=sac.EVLA; evlo=sac.EVLO;
		latmin=min([evla stadata(:).stla])-10; % Bottom Latitude
		latmax=max([evla stadata(:).stla])+10;  % Top Latitude

		longmin=min([evlo stadata(:).stlo])-10; % Left Longitude
		longmax=max([evlo stadata(:).stlo])+10; % Right Logitude
		
		if (longmax-longmin)>180
			temp=longmax;
			longmax=longmin+360+20;
			longmin=temp-20;
		end

		Pline=10; % Grid spacing for latidue lines
		Mline=10; % Grid spacing for longitude lines

		%define event location
		elat=evla; % Event Latitude
		elong=evlo; % Event Longitude

		hh=axesm('mapproj','aitoff',...
		   'maplatlim',[latmin latmax],'maplonlim',[longmin longmax],...
		   'MLineLocation',Mline,'PLineLocation',Pline,'Grid','on',...
		   'MeridianLabel','on','ParallelLabel','on');
		axis off, gridm on, framem on;
%		hh=worldmap('world');
		load coast
		plotm(lat,long);
		
		hold on

		for i=1:length(stadata)
		   [la, lo]=gcwaypts(elat,elong,stadata(i).stla,stadata(i).stlo,30);
		   hh=geoshow(la,lo,'displaytype','line','color','r');
		   set(hh,'LineWidth',1)
		end
		
		[epidist baz]=distance(mean([stadata(:).stla]),mean([stadata(:).stlo]),evla,evlo);
		stemp=sprintf('%s \n Dist=%f, Baz=%f Mag:%3f',event,epidist,baz,sac.MAG);
		title(stemp);
	% End of plot the world map

	% Inital the eventinfo structure
	stemp=sprintf('%s.win',event);
	if exist(stemp)
		fpwin=fopen(stemp,'r');
		stemp=fgetl(fpwin);
		[eventname gv1 t1 gv2 g2]=strread(stemp,'%s %f %f %f %f\n');
		eventinfo(eventnum).gv1=gv1;
		eventinfo(eventnum).gv2=gv2;
		eventinfo(eventnum).t1=t1;
		eventinfo(eventnum).t2=t2;
	else
		eventinfo(eventnum).gv1=groupv(1);
		eventinfo(eventnum).gv2=groupv(2);
		eventinfo(eventnum).t1=0;
		eventinfo(eventnum).t2=0;
	end
	eventinfo(eventnum).accept=0;

	% Interact with user
	while 1
		disp('What do you want to do?')
		disp('1: Change window beginning')
		disp('2: Change window endding')
		disp('3: Accept the result')
		disp('4: Give up this event');
		disp('5: I''m tired, output the result now');
		resp=input('','s');
		if resp=='4'
			disp('Ok, give up this event!')
			eventinfo(eventnum).accept=0;
			break;
		end
		if resp=='3'
			disp('Ok, accept this event!')
			eventinfo(eventnum).accept=1;
			stemp=sprintf('%s.win',event);
			fpout=fopen(stemp,'w');
			fprintf(fpout,'%s %f %f %f %f\n',eventinfo(eventnum).name, eventinfo(eventnum).gv1, ...
								eventinfo(eventnum).t1, eventinfo(eventnum).gv2, eventinfo(eventnum).t2);

			fclose(fpout);
			break;
		end
		if resp=='5'
			disp('Ok, Good Night!')
			sleep=1;
			break;
		end
		if resp=='1'
			figure(1)
			coor=ginput(2);
			A=[coor(1,2) 1;coor(2,2) 1];
			B=[coor(1,1);coor(2,1)];
			x=A\B;
			eventinfo(eventnum).gv1=1/x(1);
			eventinfo(eventnum).t1=x(2);
		end
		if resp=='2'
			figure(1)
			coor=ginput(2);
			A=[coor(1,2) 1;coor(2,2) 1];
			B=[coor(1,1);coor(2,1)];
			x=A\B;
			eventinfo(eventnum).gv2=1/x(1);
			eventinfo(eventnum).t2=x(2);
		end

		disp([eventinfo(eventnum).gv1 eventinfo(eventnum).t1]);
		disp([eventinfo(eventnum).gv2 eventinfo(eventnum).t2]);

		% Now, we need to flash the waveform plot by new window
		figure(1)
		clf
			hold on
			amp=(dist(end)-dist(1))/plotstanum;
			for i=1:length(dist)
				plot(stadata(i).t, stadata(i).prefdata*amp+stadata(i).dist);
				for j=1:length(periods)
					plot(stadata(i).t, stadata(i).envelop(j,:)*amp/2+stadata(i).dist,'k');
				end
			end
			gv1=eventinfo(eventnum).gv1;
			gv2=eventinfo(eventnum).gv2;
			t1=eventinfo(eventnum).t1;
			t2=eventinfo(eventnum).t2;
			plot([dist(1)/gv1+t1 dist(end)/gv1+t1],[dist(1) dist(end)],'r')
			plot([dist(1)/gv2+t2 dist(end)/gv2+t2],[dist(1) dist(end)],'r')
			for k=1:length(refgv)
				plot([dist(1)/refgv(k) dist(end)/refgv(k)],[dist(1) dist(end)],'g--')
			end
			ylim([dist(1)-ddist dist(end)+ddist])
			xlim([dist(1)/groupv(1)-500 dist(end)/groupv(2)+1000])

		% Plot filtered data
		figure(2)
		clf
		for i=1:length(periods)
			subplot(1,length(periods),i)
				hold on
				amp=(dist(end)-dist(1))/plotstanum/2;
				for j=1:length(dist)
					plot(stadata(j).t, stadata(j).fdata(i,:)*amp+stadata(j).dist);
					plot(stadata(j).t, stadata(j).envelop(i,:)*amp+stadata(j).dist,'k');
				end
				plot([dist(1)/gv1+t1 dist(end)/gv1+t1],[dist(1) dist(end)],'r')
				plot([dist(1)/gv2+t2 dist(end)/gv2+t2],[dist(1) dist(end)],'r')
				for k=1:length(refgv)
					plot([dist(1)/refgv(k) dist(end)/refgv(k)],[dist(1) dist(end)],'g--')
				end
				ylim([dist(1)-ddist dist(end)+ddist])
				xlim([dist(1)/groupv(1)-500 dist(end)/groupv(2)+1000])
				title(sprintf('%d',periods(i)));
		end

	end % End of input loop

	event=fgetl(fpev);
end % Loop of event

fclose(fpev)

if sum([eventinfo(:).accept])>0.5
	outfile=input('Please input output file name: ','s');
	while exist(outfile,'file')==2
		disp('File exist! Please change another name!')
		outfile=input('Please input output file name: ','s');
	end

	fpout=fopen(outfile,'w');
	for i=1:length(eventinfo)
		if eventinfo(i).accept==1
			fprintf(fpout,'%s %f %f %f %f\n',eventinfo(i).name, eventinfo(i).gv1, eventinfo(i).t1, eventinfo(i).gv2, eventinfo(i).t2);
		end
	end
	fclose(fpout);
end
