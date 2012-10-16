function [dz_cor,amp]=cwtgeomorlet2d_correct(dz,cwt)

	% get the real part of the cwt and normalize it.
	realcwt=real(cwt)./nanmax(nanmax(real(cwt)));
	abscwt=abs(cwt)./nanmax(nanmax(abs(cwt)));

	%calculate the average velocity of the whole map
	avgv=nanmean(nanmean(dz));

	dz0=dz-avgv;

	eng = @(a) nansum(nansum((dz0-a.*realcwt).^2.*abscwt.^2));
	%eng = @(a) nansum(nansum((dz0-a.*realcwt).^2));


	amp=fminbnd(eng,-avgv*0.2,avgv*0.2);

	dz_cor=dz-amp.*realcwt;


	


