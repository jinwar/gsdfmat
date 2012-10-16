% This code should be run after the anisotropy_averageevent.m
% written by Ge Jin
% jinwar@gmail.com
%
%



figure(9)
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
avgphv=nanmean(nanmean(newisoGV));
%seiscolormap
load seiscmap
colormap(seiscmap);
caxis([avgphv*(1-r) avgphv*(1+r)]);
u=aniso_strength.*cosd(aniso_azi)*100;
v=aniso_strength.*sind(aniso_azi)*100./cosd(40);
% Down sample the azimuth
sxnode=lalim(1):1:lalim(2);
synode=lolim(1):1:lolim(2);
[sxi,syi]=meshgrid(sxnode,synode);
su=interp2(xi,yi,u',sxi,syi);
sv=interp2(xi,yi,v',sxi,syi);
azistd=interp2(xi,yi,aniso_azi_std',sxi,syi);
strstd=interp2(xi,yi,aniso_strength_std',sxi,syi);
%         h=quiverm(sxi,syi,su,sv,'k-');
[m n]=size(sxi);
for ix=1:m
    for iy=1:n
        if azistd(ix,iy) < 40% && strstd(ix,iy)<0.01 && aniso_strength(ix,iy)>0.005
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

[mla mlo]=inputm(1);
xnode=inverse.xnode;
ynode=inverse.ynode;

mi=find(abs(xnode-mla)==min(abs(xnode-mla)));
mj=find(abs(ynode-mlo)==min(abs(ynode-mlo)));

avgV_best=isoGV_bestcor(mi,mj)
avgV=isoGV(mi,mj)

n=0;
clear phV_best azi phV dist;
for i=1:eventnum
    if ~isnan(eventdata(i).GV_bestcor(mi,mj))
        n=n+1;
        [dist(n) gcazi(n)]=distance(mla,mlo,eventdata(i).evla,eventdata(i).evlo);
        %			azi(n)=angle(eventdata(i).GVx(mi,mj) + eventdata(i).GVy(mi,mj).*sqrt(-1));
        %			azi(n)=rad2deg(azi(n));
        %			ind=find(azi<0);
        %			azi(ind)=azi(ind)+360;
        azi(n)=eventdata(i).azi(mi,mj);
        phV_best(n)=eventdata(i).GV_bestcor(mi,mj);
        phV(n)=eventdata(i).GV(mi,mj);
        besterr(n)=eventdata(i).besterr;
    end
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

%     para=fit_azi_anisotropy(azibin,aziavg_phv_best,azistd_phv_best);
para2=fit_azi_anisotropy(azi,phV_best);
para=fit_azi_anisotropy_1phi(azi,phV_best);

figure(98)
clf
hold on
plot(azi,phV,'x');
%         h1=plot(para,'k');
%         h2=plot(para2,'b');
%         legend(h1,'1phi',h2,'2phi');
% 		plot([0 360],[avgV avgV]);
% plot([0 360],[mean(phV) mean(phV)]);
% errorbar(azibin,aziavg_phv,azistd_phv,'b.')
xlim([0 360])
%         ylim([3.8 4.3])
xlabel('Azimuth')

figure(99)
clf
hold on
plot(azi,phV_best,'or');
% errorbar(azibin,aziavg_phv_best,azistd_phv_best,'r.')
h1=plot(para,'k');
% h2=plot(para2,'b');
set(h1,'linewidth',3)
% legend(h1,'1phi',h2,'2phi');
% 		plot([0 360],[avgV_best avgV_best],'r');
% plot([0 360],[mean(phV_best) mean(phV_best)]);
xlim([0 360])
avgphV_best = mean(phV_best);
ylim([avgphV_best*.9 avgphV_best*1.1]);
%         ylim([3.8 4.3])
set(gca,'fontsize',20)
xlabel('Azimuth')
ylabel('Phase Velocity (km/s)')


figure(97)
clf
hold on
errorbar(azibin,aziavg_phv,azistd_phv,'b.')
errorbar(azibin,aziavg_phv_best,azistd_phv_best,'r.')
xlim([0 360])
xlabel('Azimuth')
%	figure(99)
%	clf
%	hold on
%		plot(dist,phV,'x');
%		plot([0 180],[avgV avgV]);
%		xlabel('Distance')
