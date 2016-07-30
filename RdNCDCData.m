function [dttm,timemin,wnddatenum,wndspeed,wnddir,pres,temp]=RdNCDCData(filename,mo,yr,iminsamp)

% RdPFWndData reads one month of wind speed info from NCDC data file
% RdPFWndData(filename, mo, yr) opens a file containing airpoort weather data
% and reads in a wind speed as well as when the data was obtained

%    J. W. Naughton
%    Department of Mechanical Engineering
%    University of Wyoming 
%    University of Wyoming Aeronautical Laboratory
%    (307)766-6284
%
%      Date     Name           Modification 
%      1/14     Naughton       Original      
%
%     Function arguments
%       filename - file to read
%       mo - month to read
%       yr -year to read
%       iminsamp - the minute of the hour to read to avoid multiple
%                   readings in one hour
%
%     Function return values
%       dttm - date and time data taken
%       timmin - time in minutes from the beginning of the month
%       wnddatenum - date/time in serial date number
%       wndspeed - wndspeed in m/s
%       wnddir - wind direction in degrees, -1 indicates no direction
%       pres - pressure in Pa
%       temp - array of temperatures in degrees Celsius
 
%****************    Execution Begins     *****************************  

%clear all
%close all

%mo=3;
%yr=2009;
%filename='G:\Users\Bsync\Data\WindResource\NCDCData\Laramie2005_2013.dat';
%iminsamp=53;

% Preallocate arrays
%day=zeros(NRec,1);
%va=zeros(NRec,15);

fid=fopen(filename,'rt');

% Read header
tline = fgetl(fid);

% Skip data we dont want 
ii=0;
'Skipping to data requested'
while (1)
    tline=fgetl(fid);
    wsdatetime=tline(14:25);
    c=sscanf(wsdatetime,'%4d%2d%2d%2d%2d',5);
    month=c(2);
    year=c(1);
    if (month== mo && year == yr)
        break
    end
    ii=ii+1;
end
   

% Now read data we do want

i=1;
while month==mo
      
    % Save important variables
    
    % Date information
    dttm.year(i)=c(1);
    dttm.month(i)=c(2);
    dttm.day(i)=c(3);
    dttm.hour(i)=c(4);
    dttm.minute(i)=c(5);
    
    wnddatenum=datenum([c(1) c(2) c(3) c(4)  c(5)  0]);
       
    % Now store other data
    
    % Wind Speed
    cwndspeed=tline(31:33);
    if(~strcmp(cwndspeed,'***')) % Check to make sure there is good wind data
        % If so, read the rest of the data
        wndspeed(i)=str2double(cwndspeed);

        % Wind Direction
        cwnddir=tline(27:29);
        if(strcmp(cwnddir,'***')) % *** here means no direction
            wnddir(i)=NaN;
        else
            wnddir(i)=str2double(cwnddir);
        end

        % Pressure
        cpres=tline(107:112);
        pres(i)=str2double(cpres);

        %Temperature
        ctemp=tline(84:88);
        temp(i)=str2double(ctemp);
    end
    
    % Read the next line
    tline=fgetl(fid);
    if(isempty(tline))
       'got here'
        break
    end

    wsdatetime=tline(14:25);
    c=sscanf(wsdatetime,'%4d%2d%2d%2d%2d',5);
    month=c(2);
    year=c(1);

    
    i=i+1;
    ii=ii+1;
    [i ii]
end

fclose(fid);

% Determine the time of each point
timemin=dttm.day*1440+dttm.hour*60+dttm.minute;

% Keep only windspeeds that are returned at indicated number past the hour

indx=(dttm.minute>iminsamp-3 & dttm.minute<iminsamp+3);
wndspeed=wndspeed(indx)*0.447;  % Convert from mph to m/s
wnddir=wnddir(indx);
pres=pres(indx)*100; % Convert from millibar to Pa
temp=temp(indx);

timemin=timemin(indx);


