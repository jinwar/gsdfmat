function d=dist(lat1,lon1,lat2,lon2)
% This function is used to calculate the rough distance between two coors quickly.

  theta = lon1 - lon2; 
  d = sin(degtorad(lat1)) .* sin(degtorad(lat2)) + cos(degtorad(lat1)).*cos(degtorad(lat2)).*cos(degtorad(theta)); 
  d = acos(d); 
  d = d * 6371; 

