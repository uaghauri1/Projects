function pbr2
clc; close all; clear all;
%this is Usman Ghauri's personal contribution to the project

T0=650; %inlet temperaturein K
P0=30*101; %inlet pressure in kPa
R=8.314; %gas constant in LkPa/molk
ft=37.07;
Fp0=35.3665; %inlet flowrate of propylene in gmol/s
Fb0=45.4276; %inlet flowrate of benzene ingmol/s 
Fi=1.7015; %flowrate of inert in gmol/s
phi=0.5;
g=9.81;



%Initial Conditions [propylene;benzene;pdipb;cumene;temp;inert]
F0=[Fp0;Fb0;0;0.0016;T0;Fi]; 

w=0:100:1.2e6; % in g
[wnew,Fnew]=ode45(@MIMO_Exp,w,F0);


FP=Fnew(:,1); %Propylene Flowrates
Xp=(-FP+Fp0)/Fp0; %Conversion

FD=Fnew(:,3); %p-DIPB

FC=Fnew(:,4); %cumene

F12=[Fnew(:,1) Fnew(:,2) Fnew(:,3) Fnew(:,4)]; %all the flowrates


figure(1);
plot(w,Xp,'linewidth',2),grid
title('Conversion')
xlabel('Catalyst Weight in g'), ylabel('Conversion of Propylene')
 
figure(2);
hold on
yyaxis left
plot(wnew,F12,'linewidth',1);
title('Molar Flow rates and Temperature')
xlabel('w'); ylabel('Flowrates in mol/s'); legend('Fp','Fb','Fd','Fc');

yyaxis right
plot(wnew,Fnew(:,5));
ylabel('Temperature (K)')
hold off


 
 

function Fdot=MIMO_Exp(w,F)
%Below are the A,B,C,D values used for the Cp values. 
Ai=-4.224; Bi=0.3063; Ci=-1.586*10^-4; Di=3.215*10^-8; % i is inert
Ap=3.71; Bp=0.2345; Cp=-1.16*10^-4; Dp=2.205*10^-8;
Ab=-33.92; Bb=0.4739; Cb=-3.017*10^-4; Db=7.130*10^-8;
Ac=-41.15; Bc=0.7935; Cc=-0.0005; Dc=1*10^-7;
Ad=-63.668; Bd=1.1737; Cd=-8*10^-4; Dd=2.1*10^-7;

%Standard Heat of Formations in J/gmol
Hfi=-104700; Hfp=20430; Hfb=82980; Hfc=3940; Hfd=-76510; 

Hr1=Hfc-(Hfp+Hfb); Hr2=Hfd-(Hfc+Hfp); %Heat of reactions at standard conditions
Fi=6.383; %inert
R=8.314; %LkPa/Kmol
Ua=9.92e14; % to be used when non adiabatic
Ta=298.15; %K
P0=25*101; %kPa

 


Fp=F(1); Fb=F(2); Fd=F(3); Fc=F(4); T=F(5); Fi=F(6);

rhocat=2000; %bulk density of catalyst
% HR=U*a*(Ta-T)/rhocat;
Ft=Fp+Fb+Fc+Fi+Fd; %total flowrates in gmol/s
k1=3.5*(10^4)*exp(-24.9*4186.8/(R*T));  %L^2/mol/s/gcat
k2=2.9*(10^6)*exp(-35.08*4186.8/(R*T)); %L^2/mol/s/gcat

%heat generated
HG=(k1*(P0^2)*Fp*Fb/(R^2*T^2*Ft^2))*(-Hr1+(T-Ta)*(Ap+Ab-Ac)+(T^2-Ta^2)*((Bp+Bb-Bc)/2)+(T^3-Ta^3)*((Cp+Cb-Cc)/3)+(T^4-Ta^4)*((Dp+Db-Dc)/4))+(k2*P0^2*Fc*Fp/(R^2*T^2*Ft^2))*(-Hr2+(T-Ta)*(Ap+Ac-Ad)+(T^2-Ta^2)*((Bp+Bc-Bd)/2)+(T^3-Ta^3)*((Cp+Cc-Cd)/3)+(T^4-Ta^4)*((Dp+Dc-Dd)/4));

%denominator for dT/dW
DEN=Fp*(Ap+Bp*T+Cp*T^2+Dp*T^3)+Fb*(Ab+Bb*T+Cb*T^2+Db*T^3)+Fc*(Ac+Bc*T+Cc*T^2+Dc*T^3)+Fd*(Ad+Bd*T+Cd*T^2+Dd*T^3)+Fi*(Ai+Bi*T+Ci*T^2+Di*T^3);


f1=((P0^2)/((R^2)*(T^2)*(Ft^2)))*(-k1*Fp*Fb-k2*Fp*Fc); %Propylene

f2=(-k1*(P0^2)*Fp*Fb)/(R^2*T^2*Ft^2); %Benzene

f3=(k2*(P0^2)*Fp*Fc)/(R^2*T^2*Ft^2); %pDIPB

f4=((P0^2)/(R^2*T^2*Ft^2))*(k1*Fp*Fb-k2*Fp*Fc); %Cumene

f5=HG/DEN;  %Temperature

f6=0; % Inert



Fdot=[f1;f2;f3;f4;f5;f6];


