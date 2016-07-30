clc
clear all
format compact

%%

filename = 'Laramie2005_2015.dat';
mo = 1;
yr = 2005;
iminsamp = 53;
[dttm,timemin,wnddatenum,wndspeed,wnddir,pres,temp]=RdNCDCData(filename,mo,yr,iminsamp)
