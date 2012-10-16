function [lat,lon]=xy2latlon(x,y,lat0,lon0)

%	x=lat-lat0;
	lat=x+lat0;
%	y=(lon-lon0).*cosd(lat);
	lon=y./cosd(lat)+lon0;
end
