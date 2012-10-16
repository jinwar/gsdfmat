% This is a function to calucate the mean without counting the zeros

function wlmean=waterlevelmean(A,wl,n)
	wlind=find(abs(A(:))>wl);
	wlmean=mean(A(wlind).^(1/n)).^n;
end
