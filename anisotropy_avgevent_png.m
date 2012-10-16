% This script is used to deal with the output of program inversestadtV6, to average the tomography from
% many event and gives the isotropy and anisotropy tomography result.
%
% Written by Ge Jin
% jinwar@gmail.com
%


clear;

%eventslist='goodcsinv0';
r=0.04;
Is1phi=1;

for period=0:7
    
    filename = sprintf('averageevent_%1d',period);
    load(filename);
    
    % Calculate the anisotropy.
    disp('Start to inverse for the azimuthal anisotropy');
    smsize=3;
    tic
    isophv=zeros(Nx,Ny);
    isophv_std=zeros(Nx,Ny);
    aniso_strength=zeros(Nx,Ny);
    aniso_strength_std=zeros(Nx,Ny);
    aniso_azi=zeros(Nx,Ny);
    aniso_azi_std=zeros(Nx,Ny);
    aniso_1phi_strength=zeros(Nx,Ny);
    aniso_1phi_azi=zeros(Nx,Ny);
    
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
                        if ~isnan(eventdata(i).GV(ii,jj))
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
            %             para=fit_azi_anisotropy(azi,phV_best);
            if Is1phi
                para=fit_azi_anisotropy_1phi(azi,phV);
            else
                para=fit_azi_anisotropy(azi,phV);
            end
            parastd=confint(para);
            isophv(mi,mj)=para.a;
            isophv_std(mi,mj)=parastd(2,1)-parastd(1,1);
            aniso_strength(mi,mj)=para.d;
            
            aniso_azi(mi,mj)=para.e;
            if para.e > 180
                aniso_azi(mi,mj)=para.e-180;
            elseif para.e < 0
                aniso_azi(mi,mj)=para.e+180;
            end
            if Is1phi
                aniso_1phi_strength(mi,mj)=para.b;
                aniso_1phi_azi(mi,mj)=para.c;
                aniso_strength_std(mi,mj)=parastd(2,4)-parastd(1,4);
                aniso_azi_std(mi,mj)=parastd(2,5)-parastd(1,5);
            else
                aniso_strength_std(mi,mj)=parastd(2,2)-parastd(1,2);
                aniso_azi_std(mi,mj)=parastd(2,3)-parastd(1,3);
            end
            
        end
    end % end of looping the grids.
    toc
    %     ind=find(aniso_azi > 180);
    %     aniso_azi(ind)=aniso_azi(ind)-180;
    %     ind=find(aniso_azi < 0);
    %     aniso_azi(ind)=aniso_azi(ind)+180;
    filename = sprintf('anisotropy_%1d',period);
    save(filename);
    
    %     figure(1)
    %     clf
    %     hold on
    %     ax = worldmap(lalim, lolim);
    %     set(ax, 'Visible', 'off')
    %     states = shaperead('usastatehi', 'UseGeoCoords', true);
    %     geoshow(ax, states, 'FaceColor', [0.7 0.7 0.7])
    % %     surfacem(xi,yi,newisoGV')
    %     surfacem(xi,yi,newisoGV')
    %     geoshow(ax, states, 'FaceColor', 'none')
    %     %		plotm(stadata(:,2),stadata(:,3),'rv');
    %     %contourm(xi,yi,z,30,'k');
    %     %seiscolormap
    %     load seiscmap
    %     colormap(seiscmap);
    %     colorbar
    %     avgphv=nanmean(nanmean(newisoGV));
    %     caxis([avgphv*(1-r) avgphv*(1+r)]);
    %     filename = sprintf('isomap_grd_%3d',periods(period+1));
    %     title(filename,'Interpreter','none');
    %     filename = sprintf('avg_isomap_grd_%1d',period);
    %     print('-djpeg99',filename);
    %
    %     figure(2)
    %     clf
    %     hold on
    %     ax = worldmap(lalim, lolim);
    %     set(ax, 'Visible', 'off')
    %     states = shaperead('usastatehi', 'UseGeoCoords', true);
    %     geoshow(ax, states, 'FaceColor', [0.7 0.7 0.7])
    %     surfacem(xi,yi,newisoGV_cor')
    %     %		plotm(stadata(:,2),stadata(:,3),'rv');
    %     %contourm(xi,yi,z,30,'k');
    %     geoshow(ax, states, 'FaceColor', 'none')
    %     %seiscolormap
    %     load seiscmap
    %     colormap(seiscmap);
    %     caxis([avgphv*(1-r) avgphv*(1+r)]);
    %     colorbar
    %     filename = sprintf('isomap_cor_%3d',periods(period+1));
    %     title(filename,'Interpreter','none');
    %     filename = sprintf('avg_isomap_cor_%1d',period);
    %     print('-djpeg99',filename);
    %
    %     figure(3)
    %     clf
    %     hold on
    %     ax = worldmap(lalim, lolim);
    %     set(ax, 'Visible', 'off')
    %     states = shaperead('usastatehi', 'UseGeoCoords', true);
    %     surfacem(xi,yi,bestcount')
    %     %		plotm(stadata(:,2),stadata(:,3),'rv');
    %     %contourm(xi,yi,z,30,'k');
    %     geoshow(ax, states, 'FaceColor', 'none')
    %     %seiscolormap
    %     load seiscmap
    %     colormap(seiscmap);
    %     colorbar
    %     filename = sprintf('datadensity_%3d',periods(period+1));
    %     title(filename,'Interpreter','none');
    %     filename = sprintf('avg_datadensity_%1d',period);
    %     print('-djpeg99',filename);
    % %
    %     figure(4)
    %     clf
    %     hold on
    %     ax = worldmap(lalim, lolim);
    %     set(ax, 'Visible', 'off')
    %     states = shaperead('usastatehi', 'UseGeoCoords', true);
    %     geoshow(ax, states, 'FaceColor', [0.7 0.7 0.7])
    %     surfacem(xi,yi,isoGV_bestcor')
    %     %		plotm(stadata(:,2),stadata(:,3),'rv');
    %     %contourm(xi,yi,z,30,'k');
    %     geoshow(ax, states, 'FaceColor', 'none')
    %     %seiscolormap
    %     load seiscmap
    %     colormap(seiscmap);
    %     caxis([avgphv*(1-r) avgphv*(1+r)]);
    %     colorbar
    %     filename = sprintf('isomap_bestcor_%3d',periods(period+1));
    %     title(filename,'Interpreter','none');
    %     filename = sprintf('avg_isomap_bestcor_%1d',period);
    %     print('-djpeg99',filename);
    
    figure(5)
    clf
    hold on
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')
    surfacem(xi,yi,isophv')
    drawpng
    %		plotm(stadata(:,2),stadata(:,3),'rv');
    %seiscolormap
    for i=1:2:eventnum
        plotm(eventdata(i).stadata(:,2),eventdata(i).stadata(:,3),'rv');
    end

    load seiscmap
    colormap(seiscmap);
    avgphv=nanmean(nanmean(isophv));
    caxis([avgphv*(1-r) avgphv*(1+r)]);
    colorbar
    filename = sprintf('anisomap_bestcor_%3d',periods(period+1));
    title(filename,'Interpreter','none');
    filename = sprintf('anisophv_bestcor_%1d',period);
    print('-djpeg99',filename);
    
    % 		figure(6)
    % 		clf
    % 		hold on
    % 		ax = worldmap(lalim, lolim);
    % 		set(ax, 'Visible', 'off')
    % 		states = shaperead('usastatehi', 'UseGeoCoords', true);
    % 		geoshow(ax, states, 'FaceColor', [0.7 0.7 0.7])
    % 		surfacem(xi,yi,isophv_std')
    % 		%		plotm(stadata(:,2),stadata(:,3),'rv');
    % 		%contourm(xi,yi,z,30,'k');
    % 		geoshow(ax, states, 'FaceColor', 'none')
    % 		colorbar
    % 		filename = sprintf('anisophv_std_bestcor_%3d',periods(period+1));
    % 		title(filename,'Interpreter','none');
    % 		filename = sprintf('anisophv_std_bestcor_%1d',period);
    % 		print('-djpeg99',filename);
    
    figure(7)
    clf
    hold on
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')
    states = shaperead('usastatehi', 'UseGeoCoords', true);
    geoshow(ax, states, 'FaceColor', [0.7 0.7 0.7])
    surfacem(xi,yi,aniso_strength')
    % 		caxis([0 0.02]);
    %		plotm(stadata(:,2),stadata(:,3),'rv');
    %contourm(xi,yi,z,30,'k');
    geoshow(ax, states, 'FaceColor', 'none')
    colorbar
    filename = sprintf('aniso_strength_bestcor_%3d',periods(period+1));
    title(filename,'Interpreter','none');
    filename = sprintf('anisostr_bestcor_%1d',period);
    print('-djpeg99',filename);
    %
    % 		figure(8)
    % 		clf
    % 		hold on
    % 		ax = worldmap(lalim, lolim);
    % 		set(ax, 'Visible', 'off')
    % 		states = shaperead('usastatehi', 'UseGeoCoords', true);
    % 		geoshow(ax, states, 'FaceColor', [0.7 0.7 0.7])
    % 		surfacem(xi,yi,aniso_strength_std')
    % % 		caxis([0 0.03]);
    % 		%		plotm(stadata(:,2),stadata(:,3),'rv');
    % 		%contourm(xi,yi,z,30,'k');
    % 		geoshow(ax, states, 'FaceColor', 'none')
    % 		colorbar
    % 		filename = sprintf('aniso_strength_std_bestcor_%3d',periods(period+1));
    % 		title(filename,'Interpreter','none');
    % 		filename = sprintf('anisostr_std_bestcor_%1d',period);
    % 		print('-djpeg99',filename);
    %
    %
    figure(9)
    clf
    hold on
    ax = worldmap(lalim, lolim);
    set(ax, 'Visible', 'off')
    % 		states = shaperead('usastatehi', 'UseGeoCoords', true);
    % 		geoshow(ax, states, 'FaceColor', [0.7 0.7 0.7])

    surfacem(xi,yi,isophv')
    drawpng
    %		plotm(stadata(:,2),stadata(:,3),'rv');
    %contourm(xi,yi,z,30,'k');
    geoshow(ax, states, 'FaceColor', 'none')
    avgphv=nanmean(nanmean(isophv));
    %seiscolormap
    load seiscmap
    colormap(seiscmap);
    caxis([avgphv*(1-r) avgphv*(1+r)]);
    u=aniso_strength.*cosd(aniso_azi)*10;
    v=aniso_strength.*sind(aniso_azi)*10./cosd(8);
    %         u=smoothmap(xi,yi,u',100);
    %         v=smoothmap(xi,yi,v',100);
    % Down sample the azimuth
    sxnode=lalim(1):0.3:lalim(2);
    synode=lolim(1):0.3:lolim(2);
    [sxi,syi]=meshgrid(sxnode,synode);
    su=interp2(xi,yi,u',sxi,syi);
    sv=interp2(xi,yi,v',sxi,syi);
    azistd=interp2(xi,yi,aniso_azi_std',sxi,syi);
    strstd=interp2(xi,yi,aniso_strength_std',sxi,syi);
    %         h=quiverm(sxi,syi,su,sv,'k-');
    [m n]=size(sxi);
    for ix=1:m
        for iy=1:n
            if azistd(ix,iy) < 40 && strstd(ix,iy)<0.01
                h=plotm([sxi(ix,iy)-su(ix,iy)/2 sxi(ix,iy)+su(ix,iy)/2],...
                    [syi(ix,iy)-sv(ix,iy)/2 syi(ix,iy)+sv(ix,iy)/2],'k-');
                set(h,'linewidth',2)
                
            end
        end
    end
    colorbar
    %         set(h,'ShowArrowHead','off');
    filename = sprintf('aniso_azi_bestcor_%3d',periods(period+1));
    title(filename,'Interpreter','none');
    filename = sprintf('aniso_azi_bestcor_%1d',period);
    print('-djpeg99',filename);
    
    if Is1phi
        figure(10)
        clf
        hold on
        ax = worldmap(lalim, lolim);
        set(ax, 'Visible', 'off')
        states = shaperead('usastatehi', 'UseGeoCoords', true);
        geoshow(ax, states, 'FaceColor', [0.7 0.7 0.7])
        surfacem(xi,yi,isophv')
        drawpng
        %		plotm(stadata(:,2),stadata(:,3),'rv');
        %contourm(xi,yi,z,30,'k');
        geoshow(ax, states, 'FaceColor', 'none')
        avgphv=nanmean(nanmean(newisoGV));
        %seiscolormap
        load seiscmap
        colormap(seiscmap);
        caxis([avgphv*(1-r) avgphv*(1+r)]);
        u=aniso_1phi_strength.*cosd(aniso_1phi_azi)*10;
        v=aniso_1phi_strength.*sind(aniso_1phi_azi)*10./cosd(8);
        %         u=smoothmap(xi,yi,u',150);
        %         v=smoothmap(xi,yi,v',150);
        % Down sample the azimuth
        sxnode=lalim(1):0.3:lalim(2);
        synode=lolim(1):0.3:lolim(2);
        [sxi,syi]=meshgrid(sxnode,synode);
        su=interp2(xi,yi,u',sxi,syi);
        sv=interp2(xi,yi,v',sxi,syi);
        azistd=interp2(xi,yi,aniso_azi_std',sxi,syi);
        strstd=interp2(xi,yi,aniso_strength_std',sxi,syi);
        h=quiverm(sxi,syi,su,sv,'k-');
        [m n]=size(sxi);
        for ix=1:m
            for iy=1:n
                if azistd(ix,iy) < 40 && strstd(ix,iy)<0.01
                    h=plotm([sxi(ix,iy)-su(ix,iy)/2 sxi(ix,iy)+su(ix,iy)/2],...
                        [syi(ix,iy)-sv(ix,iy)/2 syi(ix,iy)+sv(ix,iy)/2],'k-');
                    set(h,'linewidth',3)
                end
            end
        end
        colorbar
        %         set(h,'ShowArrowHead','off');
        filename = sprintf('aniso_azi_bestcor_%3d',periods(period+1));
        title(filename,'Interpreter','none');
        filename = sprintf('aniso_1phi_azi_bestcor_%1d',period);
        print('-djpeg99',filename);
        print('-djpeg99',filename);
        
        figure(11)
        clf
        hold on
        ax = worldmap(lalim, lolim);
        set(ax, 'Visible', 'off')
        states = shaperead('usastatehi', 'UseGeoCoords', true);
        geoshow(ax, states, 'FaceColor', [0.7 0.7 0.7])
        surfacem(xi,yi,aniso_1phi_strength')
                drawpng
        % 		caxis([0 0.03]);
        %		plotm(stadata(:,2),stadata(:,3),'rv');
        %contourm(xi,yi,z,30,'k');
        geoshow(ax, states, 'FaceColor', 'none')
        colorbar
        filename = sprintf('aniso_1phi_str_bestcor_%3d',periods(period+1));
        title(filename,'Interpreter','none');
        filename = sprintf('aniso_1phi_str_bestcor_%1d',period);
        print('-djpeg99',filename);
        
    end % end of is1phi plot figure 11
    
end % end of period loop
