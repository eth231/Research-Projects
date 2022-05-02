clc
clear all
%% Residual Stress Plot for gamma-TiAl to find out the best cutting speed
b = 0.00006;
vC = 1.4;
Fc = 40;
rho = 3750;
T = 1:1:1000;
E = 144000;
alpha = -1.7E-28.*T.^3 - 5.4E-12.*T.^2 + 1.14E-8.*T + 1.09E-5;
k = 14.09;
cP = 602.38;
w = 0.019;
Pecletnum = (vC*b*rho*cP/(2*k));
c4 = 0.00527*Pecletnum^3 - 0.192*Pecletnum^2 + 2.39*Pecletnum;
if(Pecletnum > 5)
    FlashTemp = ((0.399*2*Fc*vC)/(k*w))*sqrt(k/(rho*cP*b));
elseif(Pecletnum < 5)
    FlashTemp = c4*0.159*(2*Fc/(rho*cP*w*b));
end
%% Normalized Von Mises Map for Frictional Sliding under Plane Strain Conditions

z_a = 0:0.1:4;
x_a = -5:0.25:5;
mu = Fc/75;       % Friction Coefficient
for i = 1:size(z_a,2)
    for j = 1:size(x_a,2)
m_a(i, j) = 0.5*((sqrt(1 + (z_a(1, j)^2)^2 - 2*(1 + (z_a(1,j))^2)*(x_a(1, i))^(2) + (x_a(1,i))^4 + 4*(z_a(1, j))^(2)*(x_a(1, i))^2)) + 1 + (z_a(1, j))^2 - (x_a(1, i))^2);
n_a(i,j) = 0.5*((sqrt((1 + (z_a(1, j))^2)^2 - 2*(1 + (z_a(1, j))^2)*(x_a(1, i))^(2) + (x_a(1, i))^4 + 4*(z_a(1, j))^(2)*(x_a(1, i))^2)) -1 - (z_a(1, j))^2 + (x_a(1, i))^2);
    end
end

for i = 1:size(z_a,2)
    for j = 1:size(x_a,2)
        sigmax1(i, j) = -(sqrt(m_a(i, j))*((1 + ((z_a(1, j)^2) + n_a(i, j)/(m_a(i, j) + n_a(i, j)))) - 2*z_a(1, j)));
        sigmax2(i, j) = -(mu*sign(x_a(1, i))*(sqrt(n_a(i, j)))*(2 - ((z_a(1, j)^2 - m_a(i, j))/(m_a(i, j) + n_a(i, j))) -2*x_a(1, i)));
        sigmax(i, j) = sigmax1(i, j) + sigmax2(i, j);
        sigmaz1(i, j) = -(sqrt(m_a(i,j))*(1 - ((z_a(1,j)^2) + n_a(i,j))/((m_a(i,j) + n_a(i,j)))));
        sigmaz2(i, j) = -mu*(sqrt(n_a(i,j))*((m_a(i,j) - (z_a(1, j)^2)/m_a(i,j) + n_a(i,j))));
        sigmaz(i, j) = sigmaz1(i,j) + sigmaz2(i,j);
    end
end
%% plot(x_a, sigmax);
v = 0.24;
YS = -9.4E-7.*T.^3 + 0.00059.*T.^2 - 0.724.*T + 720;
sigmae = (E.*alpha.*(T))/(1 - v);
for i = T
    if sigmae(i) > YS(i)
        sigmat(i) = sigmae(i) - YS(i);
    else 
        sigmat(i) = 0;
end
end 
title('Thermal Residual Stress');
xlabel('Temperature (C)');
ylabel('Stress (MPa)');
plot(T, YS);
hold on
plot(T, sigmae);
hold on
plot(T, sigmat);
title('Thermal Residual Stress');
xlabel('Temperature (^{\circ}C)');
ylabel('Stress (MPa)');
legend('Yield Stress', 'Equivalent Thermal Expansion Stress', 'Residual Stress');
