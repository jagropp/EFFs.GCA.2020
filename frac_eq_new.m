% Temperature dependent equilibrium fractionation factors
%

function a13eq = frac_eq_new
% DOCUMENTATION
% Output:
%         a13eq - a matrix of 7 alpha CO2-CH4 of the equilibrium
%         fractionation for a specific temperature
%
% Input: 
%         Tc - the temperature in celsius degrees
%
% This function provides a set of 13C equilibrium fractionation factors for
% each step, from data provided by Mark Iron for a temperature range of
% 0-500. This function extracts the data from the excel file and produces
% regression curves for each reaction. It also compared this data to
% results from Horita 2001 and Richet & Bottinga 1977.

fig_epsilon = 1;
fig_fit     = 1;
fig_sqdif   = 1;

global Tk
Tc = Tk - 273.15;

a = xlsread('Equilibrium Isotopic Fractionation HCTH',1);
temp = a(56:81,13);
beta = a(56:81,4:11); % beta values calculated by Mark
for n = 1:7
    alpha(:,n) = beta(:,n)./beta(:,n+1);
end
alpha(:,8) = beta(:,1)./beta(:,8);
epsA = 1000*log(alpha);
% eps = (alpha-1).*1000;
lnAlpha = a(56:81,14:21);

% When I have the final file I will produce a matrix that contains all the
% betas and do the calculations from there. 'beta' extracts the beta values
% from the spreadsheet. It has 8 columns each representing a compund in the
% order prescribed in the file 'methanogenesis compounds'.

if fig_epsilon
figure
plot(temp,epsA,'LineWidth',1.5);
% plot(1e6./((temp+273.15).^2),epsA,'LineWidth',1.5);
box off
legend('CO{_2} \rightarrow CHO-MFR','CHO-MFR \rightarrow CHO-H{_4}MPT',...
       'CHO-H{_4}MPT \rightarrow CH-H{_4}MPT','CH-H{_4}MPT \rightarrow CH{_2}-H{_4}MPT',...
       'CH_2-H{_4}MPT \rightarrow CH{_3}-H{_4}MPT','CH{_3}-H{_4}MPT \rightarrow CH{_3}-S-CoM',...
       'CH{_3}-S-CoM \rightarrow CH{_4}','CO{_2} \rightarrow CH{_4}','Location','northeast');
xlabel('T [{^o}C]','FontSize',14)
ylabel(['1000 ln\alpha_{CO_2-CH_4} [' char(8240) ']'],'FontSize',14)
end

% Create a polynomial fit of the 3 degree (Why third?)
for n = 1:8
   pi(n,:) = polyfit(temp,alpha(:,n),3);
   alphaFit(1,n) = pi(n,1).*(Tc.^3) + pi(n,2).*(Tc.^2) + ...
                   pi(n,3).*(Tc) + pi(n,4);
end
a13eq = alphaFit;
%plot(temp,alphaFit)

% Data from Horita 2001:
tempHor = [200.6 251.2 301.4 348.5 409 449 496.5 550.0 598.5];
lnAlphaHor = [34.24 29.33 25.19 21.98 19.26 16.93 15.29 13.87 12.7];
% The function to fit his data (Eq. 6 in Horita):
Tx = (200:1:600) + 273.15; % [K] I use K because of the fitted equation
fitHor = 26.70 - 49.137e3./Tx + 40.828e6./(Tx.^2) - 7.512e9./(Tx.^3);
% plot(Tx-273.15,fitHor);

% Fixed function from Richet 1977 (Eq. 7 in Horita 2001)
Tz = (1:1:1300) + 273.15; % [K] I use K because of the fitted equation
fitRich = 0.16 + 11.754e6./Tz.^2 - 2.3655e9./Tz.^3 + 0.2054e12./Tz.^4;
% plot(Tz-273.15,fitRich);

if fig_fit
figure
plot(temp,lnAlpha(:,8),'o',...
     tempHor,lnAlphaHor,'o',...
     Tz(:,1:700)-273.15,fitRich(:,1:700),...
     'MarkerSize',8,'LineWidth',1.5);
box off
xlabel('T [{^o}C]','FontSize',14)
ylabel(['1000 ln\alpha_{CO_2-CH_4} [' char(8240) ']'],'FontSize',14)
legend('This Work',...
       'Horita 2001',...
       'Richet 1977',...
       'Location','northeast');
legend('boxoff')
end

% Produce fitted values for Horita's temp points from Richet and our 3 fits.   
fitRich2 = 0.16 + 11.754e6./(tempHor+273.15).^2 - ...
           2.3655e9./(tempHor+273.15).^3 + 0.2054e12./(tempHor+273.15).^4;

% Extraction of the list of the 1000*ln(alpha) for 3 different methods of
% calculation, HCTH, M06L and PBE (respective to 1-3 in the columns of
% 'blnAlpha'.
k = [5 6 8];
for n = 1:3
    b = xlsread('Equilibrium Isotopic Fractionation',k(n));
    blnAlpha(:,n) = b(56:81,21);
end
% save('1000lnalpha.mat','temp','blnAlpha')
% Produce a fit for the 3 methods.  
for n = 1:3
   p2(n,:) = polyfit(temp,blnAlpha(:,n),3);
   blnAlphaFit(:,n) = p2(n,1).*(tempHor.^3) + p2(n,2).*(tempHor.^2) + ...
                   p2(n,3).*tempHor + p2(n,4);
end
% Sum of squared differences
sqdif = (blnAlphaFit(1:7,:) - repmat(lnAlphaHor(1:7)',1,3)).^2;
sum_sqdif = sum(sqdif);
methods = {'HCTH' 'M06L' 'PBE'};
str1 = sprintf('%s = %.3g',methods{1},sum_sqdif(1));
str2 = sprintf('%s = %.3g',methods{2},sum_sqdif(2));
str3 = sprintf('%s = %.3g',methods{3},sum_sqdif(3));

if fig_sqdif
figure
plot(tempHor(:,1:7),lnAlphaHor(:,1:7),'*',...
     tempHor(:,1:7),fitRich2(:,1:7),'*',...     
     tempHor(:,1:7),blnAlphaFit(1:7,1),'o',...
     tempHor(:,1:7),blnAlphaFit(1:7,2),'o',...
     tempHor(:,1:7),blnAlphaFit(1:7,3),'o',...
     'MarkerSize',8,'LineWidth',2);
box off
xlabel('T [{^o}C]','FontSize',14)
ylabel(['1000 ln\alpha_{CO_2-CH_4} [' char(8240) ']'],'FontSize',14)
legend('Horita 2001',...
       'Richet 1977',...
       str1,str2,str3,...
       'Location','northeast');
legend('boxoff')
end
   

   
   
   
   