function y=isonnode(aziinput, cmt)

	azirange=5;
	y=0;
	str1=cmt(1);
	str2=cmt(4);
	azi=aziinput;

	if abs(azi-str1) < azirange
		y=1;
	elseif abs(azi-str2) < azirange
		y=1;
	end

	azi=azi+180;

	if abs(azi-str1) < azirange
		y=1;
	elseif abs(azi-str2) < azirange
		y=1;
	end

	azi=azi-360;

	if abs(azi-str1) < azirange
		y=1;
	elseif abs(azi-str2) < azirange
		y=1;
	end

end
