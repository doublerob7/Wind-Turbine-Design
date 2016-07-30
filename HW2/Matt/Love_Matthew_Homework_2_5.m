% Matthew Love
% ME 4470 Wind and Tidal Energy
% Dr. Naughton
% October 5, 2015
%                       Homework 2, Problem 4

clear      % Clears variables
clc        % Clears workspace

%% Given: 
% A commercial wind turbine has the turbine power curve given in 
% the file available on the homework web site. In another file, the 
% rotation rate of the wind turbine is also given at different wind speeds.
% The wind turbine has a rotor disk diameter of 52 m.

%% Find:
% a) Determine the coefficient of performance for this wind turbine at the 
     % velocity values listed in the data file. Assume that the data 
     % provided was acquired at sea level.
     
% b) Most often, the coefficient of performance is plotted versus the tip 
     % speed ratio (blade tip velocity divided by the incoming wind 
     % velocity). Using the rotation rate information given, re-plot the 
     % coefficient of performance as a function of tip speed. Plot the Betz 
     % and Glauert limits on this second plot. A file with the Glauert 
     % limit tabulated is provided.
     
% c) Comment on the result.
%% Solution:

% Run Homework 2 Problem 3 code to obtain data and corresponding pdf's:
    Love_Matthew_Homework_2_4();
    
% Power curve data [m/s]:
wndspeed_curve = [0,3,4.05,5.35,6.3,7.14,8.05,8.89,9.66,10.5,11.7,12.89,...
              13.77,15.38,25,25.01,50];
power_curve = [0,0,38.229,92.555,148.893,213.28,301.811,402.414,503.018,...
         607.646,732.394,814.889,845.07,85,85,0,0];
     
wndspeed_rotation = [0,4,9.45,25,25.01];
rpm = [0,14,27,27,0];

lambda = [0.0734,0.0773,0.0812,0.0852,0.0891,0.0931,0.0972,0.1012,...
    0.1053,0.1095,0.1136,0.1178,0.122,0.1262,0.1305,0.1348,0.1391,...
    0.1435,0.1479,0.1523,0.1568,0.1613,0.1659,0.1704,0.175,0.1797,...
    0.1844,0.1891,0.1939,0.1987,0.2035,0.2084,0.2134,0.2184,0.2234,...
    0.2285,0.2336,0.2387,0.244,0.2492,0.2546,0.2599,0.2654,0.2708,...
    0.2764,0.282,0.2876,0.2934,0.2991,0.305,0.3109,0.3169,0.3229,...
    0.329,0.3352,0.3415,0.3478,0.3542,0.3607,0.3673,0.3739,0.3807,...
    0.3875,0.3944,0.4014,0.4085,0.4158,0.4231,0.4305,0.438,0.4457,...
    0.4534,0.4613,0.4693,0.4774,0.4857,0.4941,0.5026,0.5113,0.5202,...
    0.5292,0.5383,0.5476,0.5571,0.5668,0.5767,0.5867,0.597,0.6075,...
    0.6182,0.6291,0.6402,0.6517,0.6633,0.6753,0.6875,0.7001,0.7129,...
    0.7261,0.7396,0.7535,0.7678,0.7825,0.7976,0.8131,0.8292,0.8457,...
    0.8628,0.8804,0.8987,0.9176,0.9372,0.9575,0.9786,1.0006,1.0235,...
    1.0473,1.0723,1.0984,1.1257,1.1545,1.1847,1.2166,1.2503,1.2861,...
    1.3241,1.3646,1.4079,1.4544,1.5045,1.5588,1.6179,1.6826,1.7539,...
    1.8331,1.9217,2.022,2.1368,2.2703,2.4282,2.6193,2.8577,3.1674,...
    3.5941,4.2387,5.3922,8.5743];

GlaurtLmt = [0.0597,0.0626,0.0656,0.0686,0.0715,0.0745,0.0774,0.0804,...
    0.0834,0.0863,0.0893,0.0922,0.0952,0.0982,0.1011,0.1041,0.107,0.110,...
    0.1129,0.1159,0.1189,0.1218,0.1248,0.1277,0.1307,0.1336,0.1366,...
    0.1396,0.1425,0.1455,0.1484,0.1514,0.1544,0.1573,0.1603,0.1632,...
    0.1662,0.1692,0.1721,0.1751,0.1781,0.1811,0.184,0.187,0.190,0.193,...
    0.196,0.1989,0.2019,0.2049,0.2079,0.2109,0.2139,0.2169,0.2199,...
    0.2229,0.2259,0.229,0.232,0.235,0.238,0.2411,0.2441,0.2471,0.2502,...
    0.2532,0.2563,0.2594,0.2624,0.2655,0.2686,0.2717,0.2748,0.2779,...
    0.281,0.2841,0.2872,0.2904,0.2935,0.2966,0.2998,0.303,0.3061,0.3093,...
    0.3125,0.3157,0.319,0.3222,0.3254,0.3287,0.3319,0.3352,0.3385,...
    0.3418,0.3451,0.3485,0.3518,0.3552,0.3586,0.3619,0.3654,0.3688,...
    0.3723,0.3757,0.3792,0.3827,0.3863,0.3899,0.3934,0.3971,0.4007,...
    0.4044,0.4081,0.4118,0.4156,0.4194,0.4232,0.4271,0.431,0.435,0.439,...
    0.443,0.4471,0.4513,0.4555,0.4598,0.4641,0.4685,0.473,0.4776,0.4822,...
    0.4869,0.4918,0.4967,0.5018,0.507,0.5123,0.5178,0.5235,0.5294,...
    0.5356,0.5421,0.5489,0.5562,0.5641,0.5728,0.5831];

% Interpolate power data:
power = interp1(wndspeed_curve,power_curve,bin_centers);

% Wind density at sea level, rho:
    rho = 1.916;                                                % [kg/m^3]

% Rotor disk swept area:
    D = 52;                                                     % [m]
    A_c = pi*(D/2)^2;                                           % [m^2]

% Dynamic power, p_u:
    P = power_curve.*1000;                                      % [W]
    
% Wind speed:
    wndspeed = wndspeed_curve;

% Coefficient of performance, C_p:
    C_p = P./(1/2*rho*wndspeed.^3*A_c);
    
% Plot results
    plot(wndspeed,C_p,'o-')
    title('Coefficient of Performance vs. Windspeed')
    xlabel('Wind Speed [m/s]')
    ylabel('Coefficient of Performance')
    
% Tip speed ratio:    
    wndspeed_new = linspace(0,30,61);
    rpm_new = interp1(wndspeed_rotation,rpm,wndspeed_new);
    linear_v = rpm_new*2*pi/60*(D/2);
    
    tip_speed_ratio = linear_v./wndspeed_new;
    P_new = interp1(wndspeed_curve,power_curve,wndspeed_new);
    
    C_p_new = P_new./(1/2*rho*wndspeed_new.^3*A_c)*1000;
    
    betz_x = linspace(0,10,100);
    betz = linspace(0.6,0.6,100);
    
    figure
    plot(tip_speed_ratio,C_p_new,'o-')
    title('Coefficient of Performance vs. Tip Speed Ratio')
    ylabel('Coefficient of Performance')
    xlabel('Tip Speed Ratio')
    hold on
    plot(lambda,GlaurtLmt,'o-')
    plot(betz_x,betz,'o-')
    legend('C_p','Glaurt','Betz')
    
 
    
    
    
    
    
    