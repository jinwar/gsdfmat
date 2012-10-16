% This script is used to deal with the output of program inversestadtV7, to average the tomography from
% many event and gives the isotropy and anisotropy tomography result.
%
% Written by Ge Jin
% jinwar@gmail.com
%


clear;

%eventslist='goodcsinv0';
eventslist='goodcs';

mineventnum=20;
Isanisotropy=0;
r=0.1;

for period=1:6
    
	period
    % read in the first event and initial the average counts and matrix
    event_fp=fopen(eventslist,'r');
    event=fgetl(event_fp);
    filename = sprintf('%s_%1d.mat',event,period);
    disp(filename);
    inverse = load(filename);
    Nx=inverse.Nx;
    Ny=inverse.Ny;
    xi=inverse.xi;
    yi=inverse.yi;
    lalim=inverse.lalim;
    lolim=inverse.lolim;
    isoGV=zeros(Nx,Ny);   % Isotropic phase velocity
    isoGV_correct=zeros(Nx,Ny);
    count=zeros(Nx,Ny);   % count of event each grid
    periods=inverse.periods;
    eventnum=0;
    
    while ischar(event)
        
        filename = sprintf('%s_%1d.mat',event,period);
        if ~exist(filename,'file')
            event=fgetl(event_fp);
            continue;
        end
		disp(filename);
        inverse = load(filename);
        eventnum=eventnum+1;
        
        for i=1:Nx
            for j=1:Ny
                if ~isnan(inverse.GV(i,j))
                    count(i,j)=count(i,j)+1;
                    isoGV(i,j)=isoGV(i,j)+inverse.GV(i,j);
                    isoGV_correct(i,j)=isoGV_correct(i,j)+inverse.GV_correct(i,j);
                end
            end
        end
%        disp(event);
        eventdata(eventnum).GV=inverse.GV;
        eventdata(eventnum).GV_correct=inverse.GV_correct;
        eventdata(eventnum).ampmap=inverse.ampmap;
        eventdata(eventnum).dAmp=inverse.dAmp;
        eventdata(eventnum).amp_term=-inverse.dAmp./inverse.ampmap./(2*pi/periods(period+1)).^2;
        %		eventdata(eventnum).sm_amp_term=inverse.amp_term;
        eventdata(eventnum).ori_sm_amp_term=inverse.amp_term;
        eventdata(eventnum).sm_amp_term=smoothmap(xi,yi, eventdata(eventnum).amp_term,100);
        eventdata(eventnum).eventid=event;
        eventdata(eventnum).evla=inverse.evla;
        eventdata(eventnum).evlo=inverse.evlo;
        GVx=inverse.GVx;
        GVy=inverse.GVy;
        azi=angle(GVx + GVy.*sqrt(-1));
        ind=find(azi<0);
        azi(ind)=azi(ind)+2*pi;
        azi=rad2deg(azi);
        eventdata(eventnum).azi=azi;
        
        event=fgetl(event_fp);
    end % end of event loop
    
    for i=1:Nx
        for j=1:Ny
            if count(i,j)<mineventnum
                isoGV(i,j)=NaN;
                isoGV_correct(i,j)=NaN;
            else
                isoGV(i,j)=isoGV(i,j)./count(i,j);
                isoGV_correct(i,j)=isoGV_correct(i,j)./count(i,j);
            end
        end
    end
    
    % Find the best ampliutude correction term multifiler
    isoGV_bestcor=zeros(Nx,Ny);
    newisoGV=zeros(Nx,Ny);
    newisoGV_cor=zeros(Nx,Ny);
    bestcount=zeros(Nx,Ny);   % count of event each grid
    
    alpha=1:0.1:10;
    for n=1:eventnum
        for i=1:length(alpha)
            eventdata(n).err(i)=ampcorerr(isoGV,eventdata(n).GV,eventdata(n).sm_amp_term,alpha(i));
        end
        besterri = find(eventdata(n).err==min(eventdata(n).err));
		if isempty(besterri)
			besterri=1;
		end
        eventdata(n).GV_bestcor = ((eventdata(n).GV).^-2 + ...
            alpha(besterri).*eventdata(n).sm_amp_term').^-.5;
        eventdata(n).besterr = eventdata(n).err(besterri);
        eventdata(n).alpha = alpha(besterri);
        if eventdata(n).besterr > 0.2
            disp('exclude event!')
            disp(['event id ',num2str(n),' event name:',eventdata(n).eventid])
            eventdata(n).GV_bestcor(:)=NaN;
        end
        
        % Sum the best corrected map
        for i=1:Nx
            for j=1:Ny
                if ~isnan(eventdata(n).GV_bestcor(i,j))
                    bestcount(i,j)=bestcount(i,j)+1;
                    isoGV_bestcor(i,j)=isoGV_bestcor(i,j) + eventdata(n).GV_bestcor(i,j);
                    newisoGV(i,j)=newisoGV(i,j) + eventdata(n).GV(i,j);
                    newisoGV_cor(i,j)=newisoGV_cor(i,j) + eventdata(n).GV_correct(i,j);
                end
            end
        end
    end
    
    % normalize the best corrected map
    isoGV_bestcor = isoGV_bestcor./bestcount;
    newisoGV = newisoGV./bestcount;
    newisoGV_cor = newisoGV_cor./bestcount;
    for i=1:Nx
        for j=1:Ny
            if bestcount(i,j)<mineventnum
                isoGV_bestcor(i,j)=NaN;
                newisoGV(i,j)=NaN;
                newisoGV_cor(i,j)=NaN;
            end
        end
    end
    
    % Do it again but with error calculated by bestcor
    isoGV_bestcor0=isoGV_bestcor;
    isoGV_bestcor=zeros(Nx,Ny);
    bestcount=zeros(Nx,Ny);   % count of event each grid
    
    alpha=1:0.1:10;
	badevent=0;
    for n=1:eventnum
        for i=1:length(alpha)
            eventdata(n).err(i)=ampcorerr(isoGV_bestcor0,eventdata(n).GV,eventdata(n).sm_amp_term,alpha(i));
        end
        besterri = find(eventdata(n).err==min(eventdata(n).err));
		if isempty(besterri)
			besterri=1;
		end
        eventdata(n).GV_bestcor = ((eventdata(n).GV).^-2 + ...
            alpha(besterri).*eventdata(n).sm_amp_term').^-.5;
        eventdata(n).besterr = eventdata(n).err(besterri);
        eventdata(n).alpha = alpha(besterri);
        if eventdata(n).besterr > 0.15
			badevent=badevent+1;
            disp('exclude event!')
            disp(['event id ',num2str(n),' event name:',eventdata(n).eventid])
            eventdata(n).GV_bestcor(:)=NaN;
        end
        
        % Sum the best corrected map
        for i=1:Nx
            for j=1:Ny
                if ~isnan(eventdata(n).GV_bestcor(i,j))
                    bestcount(i,j)=bestcount(i,j)+1;
                    isoGV_bestcor(i,j)=isoGV_bestcor(i,j) + eventdata(n).GV_bestcor(i,j);
                end
            end
        end
    end
	disp(['exclude: ',num2str(badevent),' events'])
    
    % normalize the best corrected map
    for i=1:Nx
        for j=1:Ny
            if bestcount(i,j)<mineventnum
                isoGV_bestcor(i,j)=NaN;
            else
                isoGV_bestcor(i,j)=isoGV_bestcor(i,j)./bestcount(i,j);
            end
        end
    end
    
    filename = sprintf('averageevent_%1d',period);
    save(filename);

    if Isanisotropy
    % Calculate the anisotropy.
    disp('Start to inverse for the azimuthal anisotropy');
    smsize=2;
    tic
    isophv=zeros(Nx,Ny);
    isophv_std=zeros(Nx,Ny);
    aniso_strength=zeros(Nx,Ny);
    aniso_strength_std=zeros(Nx,Ny);
    aniso_azi=zeros(Nx,Ny);
    aniso_azi_std=zeros(Nx,Ny);

    for mi=1:Nx
		disp(mi/Nx);
        for mj=1:Ny
            
            avgV_best=isoGV_bestcor(mi,mj);
            avgV=isoGV(mi,mj);
            
            n=0;
            clear phV_best azi phV dist;
            for i=1:eventnum
                lowi=max(1,mi-smsize);
                upi=min(Nx,mi+smsize);
                lowj=max(1,mj-smsize);
                upj=min(Ny,mj+smsize);
                for ii=lowi:upi
                    for jj=lowj:upj
                        if ~isnan(eventdata(i).GV_bestcor(ii,jj))
                            n=n+1;
        %                     [dist(n) gcazi(n)]=distance(mla,mlo,eventdata(i).evla,eventdata(i).evlo);
                            azi(n)=eventdata(i).azi(ii,jj);
                            phV_best(n)=eventdata(i).GV_bestcor(ii,jj);
                            phV(n)=eventdata(i).GV(ii,jj);
                            besterr(n)=eventdata(i).besterr;
                        end
                    end
                end
            end
            
            if n < mineventnum*((2*smsize).^2)
                isophv(mi,mj)=NaN;
                isophv_std(mi,mj)=NaN;
                aniso_strength(mi,mj)=NaN;
                aniso_strength_std(mi,mj)=NaN;
                aniso_azi(mi,mj)=NaN;
                aniso_azi_std(mi,mj)=NaN;
                continue;
            end
            
            
            azibin=0:20:360;
            clear azistd* aziavg*
            for n=1:length(azibin)-1
                ind=find(azi>azibin(n) & azi<azibin(n+1));
                if length(ind)>1
                    aziavg_phv_best(n)=mean(phV_best(ind));
                    aziavg_phv(n)=mean(phV(ind));
                    if length(ind)>3
                        azistd_phv_best(n)=std(phV_best(ind));
                        azistd_phv(n)=std(phV(ind));
                    else
                        azistd_phv_best(n)=besterr(n);
                        azistd_phv(n)=besterr(n);
                    end
                else
                    aziavg_phv_best(n)=NaN;
                    aziavg_phv(n)=NaN;
                    azistd_phv_best(n)=NaN;
                    azistd_phv(n)=NaN;
                end
            end
            
            azibin=azibin+(azibin(2)-azibin(1))./2;
            azibin=azibin(1:end-1);
            
%             para=fit_azi_anisotropy(azibin,aziavg_phv_best,azistd_phv_best);
            para=fit_azi_anisotropy(azi,phV_best);
			parastd=confint(para);
            isophv(mi,mj)=para.a;
            isophv_std(mi,mj)=parastd(2,1)-parastd(1,1);
            aniso_strength(mi,mj)=para.d;
            aniso_strength_std(mi,mj)=parastd(2,2)-parastd(1,2);
            aniso_azi(mi,mj)=para.e;
            if para.e > 180
                aniso_azi(mi,mj)=para.e-180;
            elseif para.e < 0
                aniso_azi(mi,mj)=para.e+180;
            end
                
            aniso_azi_std(mi,mj)=parastd(2,3)-parastd(1,3);
          
            
        end
	end % end of looping the grids.
    toc
%     ind=find(aniso_azi > 180);
%     aniso_azi(ind)=aniso_azi(ind)-180;
%     ind=find(aniso_azi < 0);
%     aniso_azi(ind)=aniso_azi(ind)+180;
    filename = sprintf('anisotropy_%1d',period);
    save(filename);

    end % end of ifanisotropy

    figure(1)
    clf
    hold on
    ax = usamap(lalim, lolim);
    set(ax, 'Visible', 'off')
    states = shaperead('usastatehi', 'UseGeoCoords', true);
    geoshow(ax, states, 'FaceColor', [0.7 0.7 0.7])
%     surfacem(xi,yi,newisoGV')
    surfacem(xi,yi,newisoGV')
    geoshow(ax, states, 'FaceColor', 'none')
    %		plotm(stadata(:,2),stadata(:,3),'rv');
    %contourm(xi,yi,z,30,'k');
    %seiscolormap
    load seiscmap
    colormap(seiscmap);
    colorbar
    avgphv=nanmean(nanmean(newisoGV));
    caxis([avgphv*(1-r) avgphv*(1+r)]);
    filename = sprintf('isomap_grd_%3d',periods(period+1));
    title(filename,'Interpreter','none');
    filename = sprintf('avg_isomap_grd_%1d',period);
    print('-dpng',filename);
    
    figure(2)
    clf
    hold on
    ax = usamap(lalim, lolim);
    set(ax, 'Visible', 'off')
    states = shaperead('usastatehi', 'UseGeoCoords', true);
    geoshow(ax, states, 'FaceColor', [0.7 0.7 0.7])
    surfacem(xi,yi,newisoGV_cor')
    %		plotm(stadata(:,2),stadata(:,3),'rv');
    %contourm(xi,yi,z,30,'k');
    geoshow(ax, states, 'FaceColor', 'none')
    %seiscolormap
    load seiscmap
    colormap(seiscmap);
    caxis([avgphv*(1-r) avgphv*(1+r)]);
    colorbar
    filename = sprintf('isomap_cor_%3d',periods(period+1));
    title(filename,'Interpreter','none');
    filename = sprintf('avg_isomap_cor_%1d',period);
    print('-dpng',filename);
    
    figure(3)
    clf
    hold on
    ax = usamap(lalim, lolim);
    set(ax, 'Visible', 'off')
    states = shaperead('usastatehi', 'UseGeoCoords', true);
    surfacem(xi,yi,bestcount')
    %		plotm(stadata(:,2),stadata(:,3),'rv');
    %contourm(xi,yi,z,30,'k');
    geoshow(ax, states, 'FaceColor', 'none')
    %seiscolormap
    load seiscmap
    colormap(seiscmap);
    colorbar
    filename = sprintf('datadensity_%3d',periods(period+1));
    title(filename,'Interpreter','none');
    filename = sprintf('avg_datadensity_%1d',period);
    print('-dpng',filename);
    
    figure(4)
    clf
    hold on
    ax = usamap(lalim, lolim);
    set(ax, 'Visible', 'off')
    states = shaperead('usastatehi', 'UseGeoCoords', true);
    geoshow(ax, states, 'FaceColor', [0.7 0.7 0.7])
    surfacem(xi,yi,isoGV_bestcor')
    %		plotm(stadata(:,2),stadata(:,3),'rv');
    %con tourm(xi,yi,z,30,'k');
    geoshow(ax, states, 'FaceColor', 'none')
    %seiscolormap
    load seiscmap
    colormap(seiscmap);
    caxis([avgphv*(1-r) avgphv*(1+r)]);
    colorbar
    filename = sprintf('isomap_bestcor_%3d',periods(period+1));
    title(filename,'Interpreter','none');
    filename = sprintf('avg_isomap_bestcor_%1d',period);
    print('-dpng',filename);
    
    if Isanisotropy

		figure(5)
		clf
		hold on
		ax = usamap(lalim, lolim);
		set(ax, 'Visible', 'off')
		states = shaperead('usastatehi', 'UseGeoCoords', true);
		geoshow(ax, states, 'FaceColor', [0.7 0.7 0.7])
		surfacem(xi,yi,isophv')
		%		plotm(stadata(:,2),stadata(:,3),'rv');
		%contourm(xi,yi,z,30,'k');
		geoshow(ax, states, 'FaceColor', 'none')
		%seiscolormap
		load seiscmap
		colormap(seiscmap);
		caxis([avgphv*0.95 avgphv*1.05]);
		colorbar
		filename = sprintf('anisomap_bestcor_%3d',periods(period+1));
		title(filename,'Interpreter','none');
		filename = sprintf('anisophv_bestcor_%1d',period);
		print('-dpng',filename);
		
		figure(6)
		clf
		hold on
		ax = usamap(lalim, lolim);
		set(ax, 'Visible', 'off')
		states = shaperead('usastatehi', 'UseGeoCoords', true);
		geoshow(ax, states, 'FaceColor', [0.7 0.7 0.7])
		surfacem(xi,yi,isophv_std')
		%		plotm(stadata(:,2),stadata(:,3),'rv');
		%contourm(xi,yi,z,30,'k');
		geoshow(ax, states, 'FaceColor', 'none')
		colorbar
		filename = sprintf('anisophv_std_bestcor_%3d',periods(period+1));
		title(filename,'Interpreter','none');
		filename = sprintf('anisophv_std_bestcor_%1d',period);
		print('-dpng',filename);
		
		figure(7)
		clf
		hold on
		ax = usamap(lalim, lolim);
		set(ax, 'Visible', 'off')
		states = shaperead('usastatehi', 'UseGeoCoords', true);
		geoshow(ax, states, 'FaceColor', [0.7 0.7 0.7])
		surfacem(xi,yi,aniso_strength')
		caxis([0 0.02]);
		%		plotm(stadata(:,2),stadata(:,3),'rv');
		%contourm(xi,yi,z,30,'k');
		geoshow(ax, states, 'FaceColor', 'none')
		colorbar
		filename = sprintf('aniso_strength_bestcor_%3d',periods(period+1));
		title(filename,'Interpreter','none');
		filename = sprintf('anisostr_bestcor_%1d',period);
		print('-dpng',filename);
		
		figure(8)
		clf
		hold on
		ax = usamap(lalim, lolim);
		set(ax, 'Visible', 'off')
		states = shaperead('usastatehi', 'UseGeoCoords', true);
		geoshow(ax, states, 'FaceColor', [0.7 0.7 0.7])
		surfacem(xi,yi,aniso_strength_std')
		caxis([0 0.03]);
		%		plotm(stadata(:,2),stadata(:,3),'rv');
		%contourm(xi,yi,z,30,'k');
		geoshow(ax, states, 'FaceColor', 'none')
		colorbar
		filename = sprintf('aniso_strength_std_bestcor_%3d',periods(period+1));
		title(filename,'Interpreter','none');
		filename = sprintf('anisostr_std_bestcor_%1d',period);
		print('-dpng',filename);
		
		
		figure(9)
		clf
		hold on
		ax = usamap(lalim, lolim);
		set(ax, 'Visible', 'off')
		states = shaperead('usastatehi', 'UseGeoCoords', true);
		geoshow(ax, states, 'FaceColor', [0.7 0.7 0.7])
		surfacem(xi,yi,aniso_azi')
	%     caxis([0 0.1]);
		%		plotm(stadata(:,2),stadata(:,3),'rv');
		%contourm(xi,yi,z,30,'k');
		geoshow(ax, states, 'FaceColor', 'none')
		colorbar
		filename = sprintf('aniso_azi_bestcor_%3d',periods(period+1));
		title(filename,'Interpreter','none');
		filename = sprintf('aniso_azi_bestcor_%1d',period);
		print('-dpng',filename);
		
    end % end of isanisotropy
    
end % end of period loop
