function v=aziphv(para,theta)
	c=para(1);
	A1=para(2);
	phi1=para(3);
	A2=para(4);
	phi2=para(5);
	v=c.*(1 + A1.*cos(theta-phi1) + A2.*cos(2*(theta-phi2)));
end

