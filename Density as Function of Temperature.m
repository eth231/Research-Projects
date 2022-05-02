clc
clear all

%% Volumetric Expansion Percentage and Density Decrease in Inconel 600
Dv = [0.35, 0.75, 1.15, 1.51, 2, 2.52, 3];
T = [100, 200, 300, 400, 500, 600, 700];
dV1 = 3*(0.000013*100) + 3*(0.000013*100)^2 + (0.000013*100)^3;
dV2 = 3*(0.000026*200) + 3*(0.000026*200)^2 + (0.000026*200)^3; 
dV3 = 3*(0.000043*200) + 3*(0.000043*200)^2 + (0.000043*200)^3;
rho1i = 8340/(1 + dV1);
rho2i = 8340/(1 + dV2);
rho3i = 8340/(1 + dV3);
k1 = 14.18;
k2 = 15.98;
k3 = 17.678;
k4 = 19.23;
k5 = 20.62;
k6 = 24.35;
k7 = 25.71;
alpha1 = 3.65;
alpha2 = 3.95;
alpha3 = 4.26;
alpha4 = 4.53;
alpha5 = 4.78;
alpha6 = 5.09;
alpha7 = 5.38;
cP1 = 0.467;
cP2 = 0.489;
cP3 = 0.503;
cP4 = 0.517;
cP5 = 0.528;
cP6 = 0.589;
cP7 = 0.592;
rho1 = k1/(alpha1*cP1);
rho2 = k2/(alpha2*cP2);
rho3 = k3/(alpha3*cP3);
rho4 = k4/(alpha4*cP4);
rho5 = k5/(alpha5*cP5);
rho6 = k6/(alpha6*cP6);
rho7 = k7/(alpha7*cP7);
rho = [rho1, rho2, rho3, rho4, rho5, rho6, rho7];
Trho = [100, 200, 300, 400, 500, 600, 700];
y1 = -1.8e-7.*Trho.^2 - 0.00026.*Trho + 8.3;
y2 = 1.5e-6.*T.^2 + 0.0032.*T + 0.029;
ax = gca;
yyaxis right
plot(T, y2,'r-o')
xlabel('Temperature Degrees-Celsius')
ylabel('Volumetric Expansion %')
ax.YColor = 'black';
hold on
yyaxis left
plot(Trho, y1, 'black-*');
legend('Density', 'Volumetric Expansion')
ylabel('Density g/cm^3)')
ax.YColor = 'black';






