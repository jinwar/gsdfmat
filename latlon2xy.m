function [x,y]=latlon2xy(lat,lon,lat0,lon0)

	x=lat-lat0;
	y=(lon-lon0).*cosd(lat);

end
