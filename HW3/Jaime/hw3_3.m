% HW3 Problem 3
clc
clear all
close all
format compact

%%
%%%%%%%%%%%%%%%%%%% HW 5 %%%%%%%%%%%%%%%%%%%
nns = 1.000:.005:1.040;
Tg = [0 -.5 -1 -1.5 -2 -2.5 -2.8 -3 -2.8]*1000; %back torque N*m
Pg = [0 -100 -200 -290 -375 -455 -500 -500 -480]*1000; %back power W
%u = 4:2:12;
%%
%for i = 1:length(u) %wind speed m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r300 = [4:2:8];
c300 = [1.5 1.4 1.3];
theta300 = [26 16 10];

r250 = [10:2:14];
c250 = [1.2 1 .8];
theta250 = [5.7 3.7 2.5];

r180 = [16:2:20];
c180 = [.6 .4 .2];
theta180 = [2 1.5 1];

rMin = 3; %m
rTip = 21; %m
omega = 30*1/60*2*pi; %rad/s
N = 3; %number of blades
rho = 1.23; %kg/m^3
u = 12; %m/s
rSeg = [r300 r250 r180];
cSeg = [c300 c250 c180];
thetaSeg = [theta300 theta250 theta180];

% import data
fID180 = fopen('DU180.dat');
fID250 = fopen('DU250.dat');
fID300 = fopen('DU300.dat');

DU180 = textscan(fID180, '%f %f %f','headerlines',1);
DU250 = textscan(fID250, '%f %f %f','headerlines',1);
DU300 = textscan(fID300, '%f %f %f','headerlines',1);
DU180 = cell2mat(DU180);
DU250 = cell2mat(DU250);
DU300 = cell2mat(DU300);

alpha180 = DU180(1:22,1);
alpha250 = DU250(1:22,1);
alpha300 = DU300(3:28,1);

Cd180 = DU180(1:22,3);
Cd250 = DU250(1:22,3);
Cd300 = DU300(3:28,3);

Cl180 = DU180(1:22,2);
Cl250 = DU250(1:22,2);
Cl300 = DU300(3:28,2);

% Step 1
a = zeros(1,length(rSeg));
aPrime = zeros(1,length(rSeg));

dela = [];
delaPrime = [];

%% Iteration loop
while all(dela)>.000001 & all(delaPrime)>.0000001
    
    a1 = a;
    aPrime1 = aPrime;
    
    phi = atand(((1-a).*u)./((1+aPrime).*omega.*rSeg));
    alpha = phi-thetaSeg;
    
    a300 = alpha(1:3);
    a250 = alpha(4:6);
    a180 = alpha(7:9);
    
    Cl_int300 = interp1(alpha300,Cl300,a300,'linear','extrap');
    Cl_int250 = interp1(alpha250,Cl250,a250,'linear','extrap');
    Cl_int180 = interp1(alpha180,Cl180,a180,'linear','extrap');
    
    Cd_int300 = interp1(alpha300,Cd300,a300,'linear','extrap');
    Cd_int250 = interp1(alpha250,Cd250,a250,'linear','extrap');
    Cd_int180 = interp1(alpha180,Cd180,a180,'linear','extrap');
    
    Cl_int = [Cl_int300 Cl_int250 Cl_int180];
    Cd_int = [Cd_int300 Cd_int250 Cd_int180];
    
    Cn = Cl_int.*cosd(phi)+Cd_int.*sind(phi);
    Ct = Cl_int.*sind(phi)-Cd_int.*cosd(phi);
    
    sig = (cSeg.*N)./(2.*pi.*rSeg); 
    
    a = (1+((4.*sind(phi).^2)./(sig.*Cn))).^-1;
    aPrime = (((4.*sind(phi).*cosd(phi))./(sig.*Ct))-1).^-1;
    dela = abs(a1-a)./a;
    delaPrime = abs(aPrime1-aPrime)./aPrime;
    
    W = u.*(1-a)./sind(phi);
    L = .5*rho.*W.^(2).*cSeg.*Cl_int;
    D = .5*rho.*W.^(2).*cSeg.*Cd_int;
    Fn = Cn.*(1/2).*rho.*W.^2.*cSeg;
    Ft = Ct.*(1/2).*rho.*W.^2.*cSeg;
    
end


%% Prandtl corrections

f = (N/2) .* (rTip-rSeg) ./ ( rSeg .* sind(phi) ) ;
F = (2/pi) .* acos( exp(-f) );

aP = ( ( (4.*F.*sind(phi).^2) ./ (sig.*Cn) ) +1 ) .^(-1);
aPrimeP = ((4.*F.*sind(phi).*cosd(phi))./(sig.*Ct)-1).^(-1);

phiP = atand(((1-aP).*u)./((1+aPrimeP).*omega.*rSeg));

alphaP = phiP - thetaSeg;

% Interpolate again here with alphaP
a300p = alphaP(1:3);
a250p = alphaP(4:6);
a180p = alphaP(7:9);

Cl_int300p = interp1(alpha300,Cl300,a300p,'linear','extrap');
Cl_int250p = interp1(alpha250,Cl250,a250p,'linear','extrap');
Cl_int180p = interp1(alpha180,Cl180,a180p,'linear','extrap');

Cd_int300p = interp1(alpha300,Cd300,a300p,'linear','extrap');
Cd_int250p = interp1(alpha250,Cd250,a250p,'linear','extrap');
Cd_int180p = interp1(alpha180,Cd180,a180p,'linear','extrap');

Cl_intP = [Cl_int300p Cl_int250p Cl_int180p];
Cd_intP = [Cd_int300p Cd_int250p Cd_int180p];

Wp = u.*(1-aP)./sind(phiP);
CnP = Cl_intP.*cosd(phiP)+Cd_intP.*sind(phiP);
CtP = Cl_intP.*sind(phiP)-Cd_intP.*cosd(phiP);
FnP = CnP .* .5 .* rho .* cSeg .* Wp.^2;
FtP = CtP .* .5 .* rho .* cSeg .* Wp.^2;


%% Part C
dT=FnP*2;
T=sum(dT).*N %Thrust, N

dTx=FtP.*rSeg.*2;
Tx=sum(dTx).*N %Torque, N-m

P=Tx.*omega %Power, W

%end

%% Plotting
figure(1)
plot(rSeg,alpha,rSeg,alphaP)
title('Angle of Attack Comparison')
xlabel('Distance from hub(m)')
ylabel('Angle of Attack (deg)')
legend('BEM','Prandtl')

figure(2)
plot(rSeg,Cn,rSeg,CnP)
title('Normal Force Coefficient Comparison')
xlabel('Distance from hub (m)')
ylabel('Normal Force Coefficient')
legend('Cn','Prandtl Cn')

figure(3)
plot(rSeg,Ct,rSeg,CtP)
title('Tangential Force Coefficient Comparison')
xlabel('Distance from Hub (m)')
ylabel('Tangential Force Coefficient')
legend('Ct','Prandtl Ct')

figure(4)
plot(rSeg,Ft,rSeg,FtP)
title('Tangential Force Comparison')
xlabel('Distance from Hub (m)')
ylabel('Tangential Force (N)')
legend('Ft','Prandtl Ft')

figure(5)
plot(rSeg,Fn,rSeg,FnP)
title('Normal Force Comparison')
xlabel('Distance from Hub (m)')
ylabel('Normal Force (N)')
legend('Fn','Prandtl Fn')


