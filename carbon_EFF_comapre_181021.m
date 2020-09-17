
fig_fit = 1;
fig_sqdif = 1;

%% 26.12.17 - EFF for carbon

[bet1,bet1lab] = xlsread('C:\Jonathan\Dropbox (Weizmann Institute)\Jonathan_(with_Itay)\Paper I EFFs\alpha-beta.xlsx',3);

Tc = bet1(:,1);
Tk = Tc + 273.15;
% gas phase
betaC(1,:) = bet1(:,4); % CO2
betaC(2,:) = bet1(:,5); % formmfr
betaC(3,:) = bet1(:,6); % formh4mpt
betaC(4,:) = bet1(:,7); % menylh4mpt
betaC(5,:) = bet1(:,8); % mleneh4mpt
betaC(6,:) = bet1(:,9); % mh4mpt
betaC(7,:) = bet1(:,10); % mcom
betaC(8,:) = bet1(:,11); % CH4

% gas phase
betaC(9,:)  = bet1(:,13); % CO2
betaC(10,:) = bet1(:,14); % formmfr
betaC(11,:) = bet1(:,15); % formh4mpt
betaC(12,:) = bet1(:,16); % menylh4mpt
betaC(13,:) = bet1(:,17); % mleneh4mpt
betaC(14,:) = bet1(:,18); % mh4mpt
betaC(15,:) = bet1(:,19); % mcom
betaC(16,:) = bet1(:,20); % CH4

for i = 1:7
    alphaC(i,:)   = betaC(i,:)./betaC(i+1,:);
    alphaC(i+7,:) = betaC(i+8,:)./betaC(i+8+1,:);
end

alphaC(15,:) = betaC(1,:)./betaC(8,:); % CO2 --> CH4 gas phase

%% 
textLab1 = {'CO_2','CHO-MFR','CHO-H_4MPT','CH-H_4MPT','CH_2-H_4MPT',...
            'CH_3-H_4MPT','CH_3-S-CoM','CH_4'};
% textLab1 = {'CO$_2$','CHO-MFR','CHO-H$_4$MPT','CH-H$_4$MPT','CH$_2$-H$_4$MPT',...
%             'CH$_3$-H$_4$MPT','CH$_3$-S-CoM','CH$_4$'};        
textLab2 = {'Fmd' 'Ftr' 'Mch' 'Mtd' 'Mer' 'Mtr' 'Mcr' 'Net'};
betaLoc = betaC(:,1);
% betaLoc(5,:) = 1.193;
% betaLoc(3,:) = 1.167;
% betaLoc(6,:) = 1.147;

col = ...
    [255,255,255;
     240,240,240;
     217,217,217;
     189,189,189;
     150,150,150;
     115,115,115;
     82,82,82;
     37,37,37;
     0,0,0]./255./1.5;
sty1 = {'-','-.','-','-.','-','-.','-','-.'};
sty2 = {'-.',':','-.',':','-.',':','-.','-'};
x1 = 1e6./((Tc+273.15).^2);
x1line = 1e6./((60+273.15).^2);
x1a = 1e6./((0+273.15).^2);
x2a = 1e6./((100+273.15).^2);

figure(1)
set(1,'Units','centimeters','Position',[10 5 20 10])
pos = [9 2 3 4 5 6 7 15];
subplot(1,2,1)
for ja = 1:8
    plot(x1,betaC((pos(ja)),:),'Color',col(ja+1,:),'LineStyle',sty1{ja})
%     text(1e3./273.15+0.1,betaLoc(ja,1),textLab1{ja},'FontName','Myriad pro')
    text(x1a+0.1,betaLoc(pos(ja),1),textLab1{ja})
    hold on
end
box off
xlabel('1000/T [^oK]')
ylabel('^{13}\beta')
% xlim([1.05 1e3./(273.15.^1)+1.5])
xlim([x2a x1a+3.3])
% ylim([0.98 1.09])
grid on
grid minor

subplot(1,2,2)
line([x1line x1line],[0.99 1.025],'Color',col(3,:),'LineWidth',0.7)
text(x1line+0.15,1.023,'60^oC','FontName','Acumin pro')
hold on
% line([1.05 1e3./(273.15.^1)+0.5],[1 1],'Color',col(3,:),'LineWidth',0.7)
% hold on
pos = 8:14;
for jb = 1:7
    plot(x1,alphaC(pos(jb),:),'Color',col(jb+1,:))
    text(x1a+0.1,alphaC(pos(jb),1),textLab2{jb})
    hold on
end
box off
xlabel('1000/T [^oK]')
ylabel('^{13}\alpha')
xlim([x2a x1a+1])
% ylim([0.98 1.09])
grid on
grid minor

%% 
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





