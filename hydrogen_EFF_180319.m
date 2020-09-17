%% 19.3.18 - EFF for hydrogen
% Comparison of hydrogen equilibrium fractionation factors for the
% H2O-H2-CH4 system.

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
col = flipud(col);
textLab1 = {'CO$_2$','CHO-MFR','CHO-H$_4$MPT','CH-H$_4$MPT','CH$_2$-H$_4$MPT',...
            'CH$_3$-H$_4$MPT','CH$_3$-S-CoM','CH$_4$'};        
textLab2 = {'Fmd' 'Ftr' 'Mch' 'Mtd' 'Mer' 'Mtr' 'Mcr' 'Mtd' 'Mer' 'Mcr' 'Frh' 'Mvh/Hdr'}; 
%% EFF for individual reactions

[a1,a2,a3] = xlsread('180221_EFF_Liu 2',10);
a1(1,:) = [];
Th = a1(:,1)';
betaH(1,:) = a1(:,2); % H2O
betaH(2,:) = a1(:,3); % formmfr
betaH(3,:) = a1(:,4); % formh4mpt
betaH(4,:) = a1(:,5); % menylh4mpt
betaH(5,:) = a1(:,6); % mleneh4mpt
betaH(6,:) = a1(:,7); % mh4mpt
betaH(7,:) = a1(:,8); % mcom
betaH(8,:) = a1(:,9); % CH4
betaH(9,:) = (a1(:,11)+a1(:,12))./2; % F420
betaH(10,:) = a1(:,13); % HS-CoB

for i = 1:7
    alphaH(i,:) = betaH(i,:)./betaH(i+1,:);
end
alphaH(8,:) = betaH(9,:)./betaH(5,:); % F420H2 -> mleneh4mpt 
alphaH(9,:) = betaH(9,:)./betaH(6,:); % F420H2 -> mh4mpt
alphaH(10,:) = betaH(10,:)./betaH(8,:); % HS-CoB -> CH4
alphaH(11,:) = betaH(1,:)./betaH(9,:); % H2O -> F420H2
alphaH(12,:) = betaH(1,:)./betaH(10,:); % H2O -> HS-CoB

alphaHLoc = alphaH(:,1);
alphaHLoc(5) = 1.07;
alphaHLoc(6) = 1.06;
alphaHLoc(7) = 1.05;
alphaHLoc(8) = 0.92;
alphaHLoc(9) = 1.06;

%% EFF from Mark
figure(1)
set(1,'Units','centimeters','Position',[10 5 20 10])
subplot(1,2,1)
x1 = 1e6./((Th+273.15).^2);
x1line = 1e6./((60+273.15).^2);
x1a = 1e6./((0+273.15).^2);
x2a = 1e6./((100+273.15).^2);
line([x1line x1line],[0.85 1.25],'Color',col(5,:),'LineWidth',0.7)
text(x1line+0.15,1.06,'60^oC','FontName','Myriad pro')
hold on
for i = 2:7
    plot(x1,alphaH(i,:),'Color',col(i,:))
    text((x1a+0.05),alphaHLoc(i),textLab2{i},'interpreter','latex')
    hold on
end
box off
% xlabel('T [^oC]')
xlabel('10^6/T^2 [K]')
ylabel('^D\alpha')
% l = legend('Fmd','Ftr','Mch','Mtd','Mer',...
%            'Mtr','Mcr',...
%            'Location','best');
% set(l,'FontSize',10)
% set(l,'interpreter','latex')
% legend('boxoff')
xlim([x2a x1a+1])
ylim([(min(min(alphaH(2:7,:))))-0.01 (max(max(alphaH(2:7,:))))+0.01])
grid on
grid minor

subplot(1,2,2)
line([x1line x1line],[0.4 2.3],'Color',col(5,:),'LineWidth',0.7)
text(x1line+0.15,2.2,'60^oC','FontName','Myriad pro')
hold on
pos = [1 8 9 10 11 12];
for i = 1:length(pos)
    plot(x1,alphaH(pos(i),:),'Color',col(i,:))
    text((x1(1)+0.05),alphaHLoc(pos(i)),textLab2{pos(i)},'interpreter','latex')
    hold on
end
box off
xlabel('T [^oC]')
xlabel('10^6/T^2 [K]')
ylabel('^D\alpha')
% l = legend('Mtd','Mer','Mcr','Frh','Mvh/Hdr',...
%            'Location','best');
% set(l,'FontSize',10)
% set(l,'interpreter','latex')
% legend('boxoff')
xlim([x2a x1a+1])
ylim([0.4 2.3])
grid on
grid minor

% %% H2O(g) - CH4
% [b1,b2,b3] = xlsread('171218_EFF_comparison',1);
% 
% T1 = b1(:,6)';
% alpha_h2o_ch4(1,:) = b1(:,8); % Suess-Horibe&Craig
% T2 = b1(:,10)';
% alpha_h2o_ch4(2,:) = b1(:,12); % Bottinga
% T3 = b1(:,14)';
% alpha_h2o_ch4(3,:) = b1(:,15); % M06L TZVP W06
% 
% figure(2)
% set(2,'Units','centimeters','Position',[10 5 18 8])
% 
% subplot(1,2,1)
% plot(T2,alpha_h2o_ch4(2,:),T3,alpha_h2o_ch4(3,:))
% box off
% xlabel('T [^oC]')
% % xlabel('10^6/T [K]')
% ylabel('\alpha')
% title('H$_2$O$_{(g)}$-CH$_4$','interpreter','latex')
% l = legend('Bottinga 1969','M06L TZVP W06',...
%            'Location','best');
% set(l,'FontSize',9)
% % set(l,'interpreter','latex')
% legend('boxoff')
% xlim([0 700])
% ylim([0.8 1.4])
% % ylim([0 3.5])
% grid on
% grid minor
% 
% subplot(1,2,2)
% plot(T2,alpha_h2o_ch4(2,:),T3,alpha_h2o_ch4(3,:),T1,alpha_h2o_ch4(1,:))
% box off
% xlabel('T [^oC]')
% % xlabel('10^6/T [K]')
% ylabel('\alpha')
% title('H$_2$O$_{(g)}$-CH$_4$','interpreter','latex')
% l = legend('Bottinga 1969','M06L TZVP W06','Suess 1949- Horibe & Craig 1995',...
%            'Location','best');
% set(l,'FontSize',9)
% % set(l,'interpreter','latex')
% legend('boxoff')
% xlim([200 500])
% ylim([0.8 1.4])
% % ylim([0 3.5])
% grid on
% grid minor
% 
% %% CH4-H2
% 
% [c1,c2,c3] = xlsread('180221_EFF_Liu 2',2);
% 
% T4 = c1(:,7)';
% alpha_ch4_h2(1,:) = c1(:,10); % Horibe&Craig
% T5 = c1(:,2)';
% alpha_ch4_h2(2,:) = c1(:,5); % Bottinga
% T6 = c1(:,12)';
% alpha_ch4_h2(3,:) = c1(:,15); % M06L TZVP W06
% 
% figure(3)
% set(3,'Units','centimeters','Position',[10 5 9 8])
% 
% % subplot(1,2,1)
% plot(T4,alpha_ch4_h2(1,:))
% hold on
% plot(T5,alpha_ch4_h2(2,:))
% hold on
% plot(T6,alpha_ch4_h2(3,:))
% box off
% xlabel('T [^oC]')
% % xlabel('10^6/T [K]')
% ylabel('\alpha')
% title('CH$_4$-H$_2$','interpreter','latex')
% l = legend('Horibe & Craig 1995','Bottinga 1969','M06L TZVP W06',...
%            'Location','best');
% set(l,'FontSize',10)
% % set(l,'interpreter','latex')
% legend('boxoff')
% xlim([0 700])
% grid on
% grid minor
% 
% %% H2O(g)-H2
% 
% [d1,d2,d3] = xlsread('171218_EFF_comparison',3);
% 
% T7 = d1(:,8)';
% alpha_h2o_h2(1,:) = d1(:,11); % Suess
% T8 = d1(:,15)';
% alpha_h2o_h2(2,:) = d1(:,18); % Bottinga
% T9 = d1(:,1)';
% alpha_h2o_h2(3,:) = d1(:,4); % M06L TZVP W06
% 
% figure(4)
% set(4,'Units','centimeters','Position',[10 5 9 8])
% 
% % subplot(1,2,1)
% plot(T7,alpha_h2o_h2(1,:))
% hold on
% plot(T8,alpha_h2o_h2(2,:))
% hold on
% plot(T9,alpha_h2o_h2(3,:))
% box off
% xlabel('T [^oC]')
% % xlabel('10^6/T [K]')
% ylabel('\alpha')
% title('H$_2$O$_{(g)}$-H$_2$','interpreter','latex')
% l = legend('Suess 1949','Bottinga 1969','M06L TZVP W06',...
%            'Location','best');
% set(l,'FontSize',10)
% % set(l,'interpreter','latex')
% legend('boxoff')
% xlim([0 700])
% grid on
% grid minor





