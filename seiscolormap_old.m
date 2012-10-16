[m,n]=size(jet(128));

matcolormap=jet(128);
for i=1:m
	seiscmap(i,:)=matcolormap(m-i+1,:);
end
save seiscmap
