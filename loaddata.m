
clear
!cat 201010200409.cs | cut -c 10-200 > cstempdata
!cat 201010200409.sta | cut -c 5-100 > statempdata
isplotmap=0;

out=load('cstempdata');
sta=load('statempdata');

