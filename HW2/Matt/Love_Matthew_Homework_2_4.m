% Matthew Love
% ME 4470 Wind and Tidal Energy
% Dr. Naughton
% October 5, 2015
%                       Homework 2, Problem 4

clear      % Clears variables
clc        % Clears workspace

%% Given/Find: The wind turbine power curve for a commercial wind turbine
% is provided in tabular form on the web site. Assume that the wind turbine
% will be installed at 80 m. Using the Weibull fit to the probability
% density function (pdf) for the wind velocity for V_e different time 
% periods (annual, January, April, July, and October) that you determined 
% in the problem above, determine and tabulate the average power and energy
% produced by the wind turbine for each of these periods. For the purpose 
% of this exercise, assume that the power curve has been adjusted for the 
% density of Laramie.

%% Solution:

% Run Homework 2 Problem 3 code to obtain data and corresponding pdf's:
    Love_Matthew_Homework_2_3();
    close all
    
% Power curve data [m/s]:
wndspeed_curve = [0,3,4.05,5.35,6.3,7.14,8.05,8.89,9.66,10.5,11.7,12.89,...
                13.77,15.38,25,25.01,50];
power_curve = [0,0,38.229,92.555,148.893,213.28,301.811,402.414,503.018,...
                607.646,732.394,814.889,845.07,850,850,0,0];
            
% Interpolate power data:
power_interpolation = interp1(wndspeed_curve,power_curve,bin_centers);

power_Jan = sum(power_interpolation.*Weibull_Jan*delta_u); % Eq. 2.83
power_Apr = sum(power_interpolation.*Weibull_Apr*delta_u);
power_Jul = sum(power_interpolation.*Weibull_Jul*delta_u);
power_Oct = sum(power_interpolation.*Weibull_Oct*delta_u);
power_year = sum(power_interpolation.*Weibull_year*delta_u);

periods = {'January';'April';'July';'October';'Year'};
Power = [power_Jan;power_Apr;power_Jul;power_Oct;power_year];
power = table(Power,'RowNames',periods);
            
    
    
    
    
    
    