% Matthew Love
% ME 4470 Wind and Tidal Energy
% Dr. Naughton
% October 5, 2015
%                       Homework 2, Problem 3

clear      % Clears variables
clc        % Clears workspace
close all

%% Given/Find: Using the 10 years of wind data provided on the web site, 
% determine the probability density function (pdf) for the wind velocity
% for each month as well as for a year. Plot your results for January,
% April, July, October and for the year.

alpha = 0.19;               % roughness coefficient (given)

%% Solution:

% Fetch desired data using RdNCDCData.m
    % Function inputs:
        filename = 'Homework_1_Laramie2005_2015.dat';
        start_mo = 1;
        start_yr = 2005;
        end_mo = 12;
        end_yr = 2005;
        iminsamp = 55;
    % Read desired data from data file:
    % Define starting point for data acquisition:
        mo = start_mo;
        yr = start_yr;
    
    % Collect and organize wind speed data
        % Love_Matthew_Extract_Data()
    % Load variables to workspace from previous session:
        load('Wind_Speed_Data.mat')

%% a) Create pdf's        
% Create histogram for each month of the year
    bin_centers = linspace(0,36,19);
    delta_u = bin_centers(2)-bin_centers(1);
    hist_Jan = hist(wndspeed_80m_Jan,bin_centers);
    hist_Apr = hist(wndspeed_80m_Apr,bin_centers);
    hist_Jul = hist(wndspeed_80m_Jul,bin_centers);
    hist_Oct = hist(wndspeed_80m_Oct,bin_centers);
    hist_year = hist(wndspeed_80m_10_yr,bin_centers);
        
% Calculate Probability Density Functions (pdf's):
    pdf_Jan = hist_Jan/(length(wndspeed_80m_Jan)*(delta_u));
    pdf_Apr = hist_Apr/(length(wndspeed_80m_Apr)*(delta_u));
    pdf_Jul = hist_Jul/(length(wndspeed_80m_Jul)*(delta_u));
    pdf_Oct = hist_Oct/(length(wndspeed_80m_Oct)*(delta_u));
    pdf_year = hist_year/(length(wndspeed_80m_10_yr)*(delta_u));
    
% Plot histogram of wind speed at 80 m during each month:
    plot(bin_centers,pdf_Jan,'o')
    title('Probability Density Function of Wind Speed Data for January')
    xlabel('Wind Speed [m/s]')
    ylabel('Probability Density Function')
    figure
    plot(bin_centers,pdf_Apr,'o')
    title('Probability Density Function of Wind Speed Data for April')
    xlabel('Wind Speed [m/s]')
    ylabel('Probability Density Function')
    figure
    plot(bin_centers,pdf_Jul,'o')
    title('Probability Density Function of Wind Speed Data for July')
    xlabel('Wind Speed [m/s]')
    ylabel('Probability Density Function')
    figure
    plot(bin_centers,pdf_Oct,'o')
    title('Probability Density Function of Wind Speed Data for October')
    xlabel('Wind Speed [m/s]')
    ylabel('Probability Density Function')
    figure
    plot(bin_centers,pdf_year,'o')
    title('Probability Density Function of Wind Speed Data for 10 Years')
    xlabel('Wind Speed [m/s]')
    ylabel('Probability Density Function')
    figure

%% b) Determine Mean Velocity and Variance
% Mean velocities
    U_Jan = sum(bin_centers.*pdf_Jan*delta_u);
    U_Apr = sum(bin_centers.*pdf_Apr*delta_u);
    U_Jul = sum(bin_centers.*pdf_Jul*delta_u);
    U_Oct = sum(bin_centers.*pdf_Oct*delta_u);
    U_year = sum(bin_centers.*pdf_year*delta_u);
    U = [U_Jan,U_Apr,U_Jul,U_Oct,U_year];
% Variances
    V_Jan = sum((bin_centers-U_Jan).^2.*pdf_Jan*delta_u);
    V_Apr = sum((bin_centers-U_Apr).^2.*pdf_Apr*delta_u);
    V_Jul = sum((bin_centers-U_Jul).^2.*pdf_Jul*delta_u);
    V_Oct = sum((bin_centers-U_Oct).^2.*pdf_Oct*delta_u);
    V_year = sum((bin_centers-U_year).^2.*pdf_year*delta_u);
    V = [V_Jan,V_Apr,V_Jul,V_Oct,V_year];
    
%% c) Weibull and Rayleigh distributions:
    Rayleigh_Jan = pi/2.*(bin_centers/U_Jan^2).*exp(-pi/4.*(bin_centers/U_Jan).^2);
    Rayleigh_Apr = pi/2.*(bin_centers/U_Apr^2).*exp(-pi/4.*(bin_centers/U_Apr).^2);
    Rayleigh_Jul = pi/2.*(bin_centers/U_Jul^2).*exp(-pi/4.*(bin_centers/U_Jul).^2);
    Rayleigh_Oct = pi/2.*(bin_centers/U_Oct^2).*exp(-pi/4.*(bin_centers/U_Oct).^2);
    Rayleigh_year = pi/2.*(bin_centers/U_year^2).*exp(-pi/4.*(bin_centers/U_year).^2);
    
% Shape factor, k:
    sigma = [sqrt(V_Jan),sqrt(V_Apr),sqrt(V_Jul),sqrt(V_Oct),sqrt(V_year)];
    k = (sigma./U).^(-1.086);
% Scale, c:
    c = U.*(0.568+0.433./k).^(-1./k);
    Weibull_Jan = (k(1)/c(1))*(bin_centers/c(1)).^(k(1)-1).*exp(-(bin_centers/c(1)).^k(1));
    Weibull_Apr = (k(1)/c(2))*(bin_centers/c(2)).^(k(2)-1).*exp(-(bin_centers/c(2)).^k(2));
    Weibull_Jul = (k(1)/c(3))*(bin_centers/c(3)).^(k(3)-1).*exp(-(bin_centers/c(3)).^k(3));
    Weibull_Oct = (k(1)/c(4))*(bin_centers/c(4)).^(k(4)-1).*exp(-(bin_centers/c(4)).^k(4));
    Weibull_year = (k(1)/c(5))*(bin_centers/c(5)).^(k(5)-1).*exp(-(bin_centers/c(5)).^k(5));
    
% Overlay Rayleigh and Weibull distributions over pdf's:
    plot(bin_centers,pdf_Jan,'o',bin_centers,Rayleigh_Jan,bin_centers,Weibull_Jan)
    title('PDF, Rayleigh and Weibull distributions: January')
    xlabel('Wind Speed [m/s]')
    ylabel('Probability Density Function')
    legend('PDF','Rayleigh','Weibull')
    figure
    
    plot(bin_centers,pdf_Apr,'o',bin_centers,Rayleigh_Apr,bin_centers,Weibull_Apr)
    title('PDF, Rayleigh and Weibull distributions: April')
    xlabel('Wind Speed [m/s]')
    ylabel('Probability Density Function')
    legend('PDF','Rayleigh','Weibull')
    figure
    
    plot(bin_centers,pdf_Jul,'o',bin_centers,Rayleigh_Jul,bin_centers,Weibull_Jul)
    title('PDF, Rayleigh and Weibull distributions: July')
    xlabel('Wind Speed [m/s]')
    ylabel('Probability Density Function')
    legend('PDF','Rayleigh','Weibull')
    figure
    
    plot(bin_centers,pdf_Oct,'o',bin_centers,Rayleigh_Oct,bin_centers,Weibull_Oct)
    title('PDF, Rayleigh and Weibull distributions: October')
    xlabel('Wind Speed [m/s]')
    ylabel('Probability Density Function')
    legend('PDF','Rayleigh','Weibull')
    figure
    
    plot(bin_centers,pdf_year,'o',bin_centers,Rayleigh_year,bin_centers,Weibull_year)
    title('PDF, Rayleigh and Weibull distributions: 10 Years')
    xlabel('Wind Speed [m/s]')
    ylabel('Probability Density Function')
    legend('PDF','Rayleigh','Weibull')