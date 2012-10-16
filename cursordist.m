function dist=cursordist()

[lat lon]=inputm(2)
dist=distance(lat(1),lon(1),lat(2),lon(2))*111;