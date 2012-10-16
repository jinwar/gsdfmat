clear;
for i=1:100
	x(i)=(rand-0.5)*100;
end

v=3;
tdt=x./v;
tv=2.7:0.01:3.3;
err=zeros(size(tv));

for i=1:length(tv)
	err(i)=sum(abs(x./tv(i)-x./v));
end

plot(err)


