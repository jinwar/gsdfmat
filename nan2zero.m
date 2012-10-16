function y=nan2zero(x)
% function to change NaN value in the matrix into zero.
% ge jin, ge.jin@ldeo

[m n]=size(x);
y=x;
for i=1:m
	for j=1:n
		if isnan(x(i,j))
			y(i,j)=0;
		end
	end
end

