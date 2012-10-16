function cmt = geteventcmt(eventname,evla,evlo)

	searchrange=2;

	data=sscanf(eventname,'%4d%2d%2d%2d%2d%2d');
	eyear=data(1);
	emonth=data(2);
	eday=data(3);
	ehour=data(4);
	emin=data(5);


	url = ['http://www.globalcmt.org/cgi-bin/globalcmt-cgi-bin/CMT4/form?itype=ymd&yr=' ...
	num2str(eyear), ...
	'&mo=' num2str(emonth) , ...
	'&day=' num2str(eday) , ...
	'&otype=ymd&oyr=' num2str(eyear), ...
	'&omo=' num2str(emonth), ...
	'&oday=' num2str(eday), ...
	'&jyr=2000&jday=1&ojyr=2010&ojday=1&nday=1', ...	% Not used, I think
	'&lmw=' '0', ...
	'&umw=' '10', ...
	'&lms=' '0', ...
	'&ums=' '10', ...
	'&lmb=' '0', ...
	'&umb=' '10', ...
	'&llat=' sprintf('%.8g', -90), ...
	'&ulat=' sprintf('%.8g', 90), ...
	'&llon=' sprintf('%.8g', -180), ...
	'&ulon=' sprintf('%.8g', 180), ...
	'&lhd=' '0', ...
	'&uhd=' '1000', ...
	'&lts=-9999&uts=9999&lpe1=0&upe1=90&lpe2=0&upe2=90&list=2']

	dest_fiche = 'lixogrr.html';
	system(['wget "' url '" -q --tries=2 --connect-timeout=5 -O ' dest_fiche]);

	finfo = dir(dest_fiche);
	if (finfo.bytes < 100)					% Delete the file anyway because it exists but is empty
		warndlg('Failed to download the CMT file.','Warning')
		builtin('delete',dest_fiche);
		set(handles.figure1,'pointer','arrow'),		return
	end

	fid = fopen(dest_fiche,'r');
	for (k = 1:24)					% Jump first 24 html related lines
		lix = fgets(fid);
	end

	todos = fread(fid,'*char');		fclose(fid);
	
	if (todos(1) == '<')
		warndlg('Could not find any event with this parameters search.','Warning')
		builtin('delete',dest_fiche);
		set(handles.figure1,'pointer','arrow'),		return
	end

	eols = find(todos == char(10));				% Find the new line breaks
	todos(eols(numel(eols)-9):end) = [];		% Remove last 9 lines that have no data

	todos';

	if (url(end) == '3')
		[lix lat lon depth lix lix mag lix xm1 xm2 xm3 xm4 xm5 xm6 ] = ...
		strread(todos,'%s %f %f %f %f %f %f %f %f %f %f %f %f %f');
		cmt=[xm1 xm2 xm3 xm4 xm5 xm6];
	else if (url(end) == '2')
		[lix lix strike dip rake str2 dip2 rake2 sc iexp name] = ...
		strread(todos,'%f %f %f %f %f %f %f %f %f %f %s');
	end
	cmt=[0 0 0 0 0 0];
	for i=1:length(name)
		temp=name(i);
		tempdata=sscanf(char(temp),'%4d%2d%2d%2d%2d%2d');
		if tempdata==data
			cmt=[strike(i) dip(i) rake(i) str2(i) dip2(i) rake2(i)];
		end
	end
	
end

