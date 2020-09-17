
fig_fit = 1;
fig_sqdif = 1;

carbon_EFF_171226

% Data from Horita 2001:
tempHor = [200.6 251.2 301.4 348.5 409 449 496.5 550.0 598.5];
epsHor = [34.24 29.33 25.19 21.98 19.26 16.93 15.29 13.87 12.7];
% The function to fit his data (Eq. 6 in Horita):
Tk1 = (200:1:600) + 273.15; % [K]
fitHor = 26.70 - 49.137e3./Tk1 + 40.828e6./(Tk1.^2) - 7.512e9./(Tk1.^3);
% plot(Tx-273.15,fitHor);

% Fixed function from Richet 1977 (Eq. 7 in Horita 2001)
Tk2 = (0:1:1573) + 273.15; % [K] I use K because of the fitted equation
fitRich = 0.16 + 11.754e6./Tk2.^2 - 2.3655e9./Tk2.^3 + 0.2054e12./Tk2.^4;
% plot(Tz-273.15,fitRich);

if fig_fit
figure
plot(tempHor,epsHor,'o',...
     Tk2-273.15,fitRich,...
     T,1000*log(alphaC(8,:)),...
     'MarkerSize',8,'LineWidth',1.5);
box off
xlabel('T [{^o}C]','FontSize',14)
ylabel(['1000 ln\alpha_{CO_2-CH_4} [' char(8240) ']'],'FontSize',14)
legend('Horita 2001','Richet 1977','M06L TZVP W06',...
       'Location','northeast');
legend('boxoff')
xlim([0 700])
end

% Produce fitted values for Horita's temp points from Richet and our 3 fits.   
fitRich2 = 0.16 + 11.754e6./(tempHor+273.15).^2 - ...
           2.3655e9./(tempHor+273.15).^3 + 0.2054e12./(tempHor+273.15).^4;

% Produce a fit for the 3 methods.  
p2 = polyfit(T,1000.*log(alphaC(8,:)),3);
lnAlphaFit = p2(1).*(tempHor.^3) + p2(2).*(tempHor.^2) + ...
    p2(3).*tempHor + p2(4);
% Sum of squared differences
sqdif = (lnAlphaFit - epsHor).^2;
sum_sqdif = sum(sqdif);
methods = {'M06L TZVP W06'};
str1 = sprintf('%s = %.3g',methods{1},sum_sqdif(1));

if fig_sqdif
figure
plot(tempHor,epsHor,'*',...
     tempHor,fitRich2,'*',...     
     tempHor,lnAlphaFit,'o',...
     'MarkerSize',8,'LineWidth',2);
box off
xlabel('T [{^o}C]','FontSize',14)
ylabel(['1000 ln\alpha_{CO_2-CH_4} [' char(8240) ']'],'FontSize',14)
legend('Horita 2001',...
       'Richet 1977',...
       str1,...
       'Location','northeast');
legend('boxoff')
xlim([min(tempHor)-50 max(tempHor)+50])
ylim([10 38])
end





