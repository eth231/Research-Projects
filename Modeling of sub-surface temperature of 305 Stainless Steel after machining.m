clc
clear all

%% Peclet Number and Flash Temperature for 305 SS

b = 40E-6;                                                   % half width of contact (m)
Fc = 40;                                                     % Cutting force in Newtons
Vc = 3*1.7;                                                      % Cutting Speed in (m/s)
w = 0.005;                                                    % Width in (m)
cPi = 500;                                                   % Specific Heat (J/(kg*C))
ki = 16.2;                                                    % Thermal Conductivity (W/(m-k))
rhoi = 7990;                                                 % density (kg/m^3)
Pe = Vc*b*rhoi*cPi/(2*ki);                                   % Peclet's Number 
TflashHigh = 0.399*2*Fc*Vc/(ki*w)*sqrt(ki/(rhoi*cPi*Vc*b));  % High Flash temperature (
C4 = 0.00527*Pe^3-0.192*Pe^2+2.39*Pe;                        % Function of Peclet's Number
TflashLow = C4*0.159*2*Fc/(rhoi*cPi*w*b);                    % Low Flash Temperature

%% If statement to determine whether Peclet's Number is small or large
if (Pe>5) 
    T = TflashHigh;
elseif (Pe<5)
    T = TflashLow;
end

% Computing Surface Temperature
x_b = -1.5:0.01:7;
Pplus = (Pe*(x_b+1))/2;
Pminus = (Pe*(x_b-1))/2;

for i=1:size(x_b,2)
    if x_b(i) <= -1 
        Tss(i) = ((x_b(i) + 1)*exp(Pplus(i))*(besselk(0,-1*Pplus(i)) - (besselk(1,-1*Pplus(i)))) + ((1 - x_b(i))*exp(Pminus(i))*(besselk(0,-1*Pminus(i)) - besselk(1,-1*Pminus(i)))));
    elseif x_b(i) > -1  && x_b(i) < 1
        Tss(i)= (((x_b(i) + 1)*exp(Pplus(i))*(besselk(0,Pplus(i)) + (besselk(1,Pplus(i))))) + ((1 - x_b(i))*exp(Pminus(i))*(besselk(0,-1*Pminus(i)) - besselk(1,-1*Pminus(i)))));
    else 
        Tss(i) = ((x_b(i) + 1)*exp(Pplus(i))*(besselk(0, Pplus(i)) + (besselk(1,Pplus(i)))) + ((1 - x_b(i))*exp(Pminus(i))*(besselk(0,Pminus(i)) + besselk(1,Pminus(i)))));
    end
end 
SurfaceTemp = ((Tss/max(Tss))*T) + 20;
grid on
title('Surface Temperature Distribution')
xlabel('x/b')
plot(x_b, SurfaceTemp);
%% Sub-Surface Temperature
Tbulk = 20;
alpha = ki/(rhoi*cPi);
t = 2*b/Vc;
z_b = 0:0.01:7;
for i= 1:size(x_b,2)
    for j = 1:size(z_b,2)      
        Tsub(i,j) = 1 - erf(z_b(1,j)*b/(sqrt((4*alpha*t)))).*SurfaceTemp(1,i) + Tbulk;
            
    end
end
plot(x_b, Tsub);
title('Thermal field of 305 SS');
xlabel('x_b')
ylabel('T');