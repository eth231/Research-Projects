clc
clear all

%Temperature Conversion Data

Temperatures = [25, 50, 100, 110, 150, 156, 181, 195, 200];
Counts = [760, 1200, 2500, 3110, 5800, 6770, 10400, 12660, 13400];
counts = 760:10:13400;
spl = spline(Counts, Temperatures);
spl2 = spl;
spl2.coefs(:,end) = spl2.coefs(:,end) + 1;
plot(Counts, Temperatures, 'o', counts, fnval(spl, counts) + 1, '-');
title('Temperature Conversion');
xlabel('Counts');
ylabel('Temperature');
hold on
plot(counts, fnval(spl2,counts), 'r--');






