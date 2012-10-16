function nanid = mapnanid(xi,yi,map,D)
	testmap = map;
	nanid = find(isnan(testmap(:)));
	vavidid = find(~isnan(testmap(:)));
	testmap(nanid)=0;
	testmap(vavidid)=1;
	sm_testmap = smoothmap_old(xi,yi,testmap,D);
	sm_testmap(vavidid)=1;
	nanid = find(sm_testmap(:)<0.5);
end
