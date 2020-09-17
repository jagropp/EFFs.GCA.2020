%% 26.12.17 - EFF for carbon

%% EFF for individual reactions

% [a1,a2,a3] = xlsread('171218_EFF_go2a',1);

load carbon_171218

Tc = a1(:,3)';
betaC(1,:) = a1(:,4); % CO2
betaC(2,:) = a1(:,5); % formmfr
betaC(3,:) = a1(:,6); % formh4mpt
betaC(4,:) = a1(:,7); % menylh4mpt
betaC(5,:) = a1(:,8); % mleneh4mpt
betaC(6,:) = a1(:,9); % mh4mpt
betaC(7,:) = a1(:,10); % mcom
betaC(8,:) = a1(:,11); % CH4

for i = 1:7
    alphaC(i,:) = betaC(i,:)./betaC(i+1,:);
end
alphaC(8,:) = betaC(1,:)./betaC(8,:); % CO2 --> CH4 

% textLab1 = {'CO_2','CHO-MFR','CHO-H_4MPT','CH-H_4MPT','CH_2-H_4MPT',...
%             'CH_3-H_4MPT','CH_3-S-CoM','CH_4'};
textLab1 = {'CO$_2$','CHO-MFR','CHO-H$_4$MPT','CH-H$_4$MPT','CH$_2$-H$_4$MPT',...
            'CH$_3$-H$_4$MPT','CH$_3$-S-CoM','CH$_4$'};        
textLab2 = {'Fmd' 'Ftr' 'Mch' 'Mtd' 'Mer' 'Mtr' 'Mcr' 'Net'};
betaLoc = betaC(:,1);
betaLoc(5,:) = 1.193;
betaLoc(3,:) = 1.167;
betaLoc(6,:) = 1.147;

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

subplot(1,2,1)
for ja = 1:8
    plot(x1,betaC((ja),:),'Color',col(ja+1,:),'LineStyle',sty1{ja})
%     text(1e3./273.15+0.1,betaLoc(ja,1),textLab1{ja},'FontName','Myriad pro')
    text(x1a+0.1,betaLoc(ja,1),textLab1{ja},'interpreter','latex')
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
line([x1line x1line],[0.98 1.1],'Color',col(3,:),'LineWidth',0.7)
text(x1line+0.15,1.085,'60^oC','FontName','Myriad pro')
hold on
% line([1.05 1e3./(273.15.^1)+0.5],[1 1],'Color',col(3,:),'LineWidth',0.7)
% hold on
for jb = 1:8
    plot(x1,alphaC((jb),:),'Color',col(jb+1,:))
%     text(1e3./273.15+0.1,alpha(jb,1),textLab2{jb},'FontName','Myriad pro')
    text(x1a+0.1,alphaC(jb,1),textLab2{jb},'interpreter','latex')
    hold on
end
box off
xlabel('1000/T [^oK]')
ylabel('^{13}\alpha')
% l = legend('Fmd','Ftr','Mch','Mtd','Mer',...
%            'Mtr','Mcr','Net',...
%            'Location','best');
% set(l,'FontSize',10)
% set(l,'interpreter','latex')
% legend('boxoff')
% xlim([1.05 1e3./(273.15.^1)+0.5])
xlim([x2a x1a+1])
ylim([0.98 1.09])
grid on
grid minor





