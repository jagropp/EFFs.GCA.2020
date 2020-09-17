clear

addpath '/Users/jonathag/Dropbox (Weizmann Institute)/Jonathan_(with_Itay)/Code/Cycling/'
extract_isotope_data
%%
figure(1)
set(1,'Units','centimeters','Position',[15 4 24 10])

cond = indLab == 1;

subplot(1,2,1)
data1 = datTable.ln_CO2_CH4(cond)'-CHdat_kln13aeq_theo(cond);
histogram(data1,20)
% title(sprintf('n = %d',sum(~isnan(data1))))
set(gca,'FontSize',14)
xlabel(['\Delta1000ln^{13}\alpha_{CO_2' char(8211) 'CH_4}'])
ylabel('n','FontAngle','italic')

subplot(1,2,2)
data2 = datTable.ln_CH4_H2O(cond)'-CHdat_kln2aeq_theo(cond);
histogram(data2,20)
% title(sprintf('n = %d',sum(~isnan(data2))))
set(gca,'FontSize',14)
xlabel(['\Delta1000ln^2\alpha_{CH_4' char(8211) 'H_2O}'])

local_path = '/Users/jonathag/Dropbox (Weizmann Institute)/Apps/Overleaf/EFF paper/figures/';
print('-f1',[local_path,'Fig_S2'],'-dpng','-r400');
print('-f1',[local_path,'Fig_S2'],'-depsc');

%%
clear
extract_isotope_data


cond2 = indLab == 4 | indLab == 5;
% figure

% subplot(1,2,1)
% data1 = datTable.ln_CO2_CH4(cond2)'-CHdat_kln13aeq_theo(cond2);
% histogram(data1,25)
% title(sprintf('n = %d',sum(~isnan(data1))))
% 
% subplot(1,2,2)
% data2 = datTable.ln_CH4_H2O(cond2)'-CHdat_kln2aeq_theo(cond2);
% histogram(data2,45)
% title(sprintf('n = %d',sum(~isnan(data2))))
%%
Tc = -10:1:100;
Tk = Tc + 273.15;
% Fit from M06-L calculations:
yM06L = -5.3937.*1e12./Tk.^4 + 53.2609.*1e9./Tk.^3 - ...
                  204.8655.*1e6./Tk.^2 + 318.7155.*1e3./Tk - 286.1796;
% Fit from HCTH calculations:
yHCTH = -5.3788.*1e12./Tk.^4 + 53.1224.*1e9./Tk.^3 - ...
                  204.2596.*1e6./Tk.^2 + 308.4225.*1e3./Tk - 283.6737;

Bottinga69aklna   = 25e6./Tk.^2 + 346e3./Tk - 223; % CH4-H2, 0 - 700 oc
Bottinga69bklna   = -7.69e6./Tk.^2 + 6.1e3./Tk + 88.4; % H2O(g)-CH4, 10-250 oc
Suess49klna       = 1000.*log(0.9714 + 225230./Tk.^2); % H2O(g)-H2(g), Eq. from Horibe 95
Richet77aklna     = 1000.*log(0.9775 + 162737./Tk.^2 + 2.9083e9./Tk.^4); % CH4-H2 (reported in Horibe95); Harmonic
Richet77bklna     = 1000.*log(0.9628 + 148618./Tk.^2 + 2.2442e9./Tk.^4); % CH4-H2 (reported in Horibe95); Anharmonic
Richet77cklna     = 4.4480e9./Tk.^3 - 46.1617e6./Tk.^2 + 131.9533e3./Tk - 8.9370; % H2O(g)-CH4; Harmonic
Horibe95Calcklna  = 1000.*log(0.8994 + 183540./Tk.^2); % CH4-H2, extrapolation of experimental data
Cerrai54Calcklna  = 1000.*log(1.045 + 211286./Tk.^2); % H2O(g)-H2(g), Eq. from Horibe 95
Rolston76Calcklna = 27.87e6./Tk.^2 + 368.9e3./Tk - 214.3; % H2O(l)-H2, calibrated from exp
Bardo76Calcklna   = 1000.*(-0.17739 + 1.12486.*(300./Tk) + 0.61704.*(300./Tk).^2 - ...
                        0.36997.*(300./Tk).^3 + 0.07978.*(300./Tk).^4); % H2O(g)-H2(g), Bardo and Wolfsberg
Horita94_h2o_l_g_klna = -0.353e18./Tk.^6 + 42.17e12./Tk.^4 - ...
    309.4e9./Tk.^3 + 963.7e6./Tk.^2 - 1399e3./Tk + 766.2;

% Extrapolation to Horibe and Craig 95
yHC95_1 = 1000.*log(0.8994 + 183540./Tk.^2)-Rolston76Calcklna; % CH4-H2, extrapolation of experimental data
yHC95_2 = 1000.*log(0.8994 + 183540./Tk.^2)-Suess49klna-Horita94_h2o_l_g_klna; % CH4-H2, extrapolation of experimental data

% Fit from M06-L calculations:
yM06L_data = -5.3937.*1e12./CHdat_Tk(cond2).^4 + 53.2609.*1e9./CHdat_Tk(cond2).^3 - ...
                  204.8655.*1e6./CHdat_Tk(cond2).^2 + 318.7155.*1e3./CHdat_Tk(cond2) - 286.1796;
% Fit from HCTH calculations:
yHCTH_data = -5.3788.*1e12./CHdat_Tk(cond2).^4 + 53.1224.*1e9./CHdat_Tk(cond2).^3 - ...
                  204.2596.*1e6./CHdat_Tk(cond2).^2 + 308.4225.*1e3./CHdat_Tk(cond2) - 283.6737;

figure(2)              
set(2,'Units','centimeters','position',[2 10 22 12.5])
ax(1) = axes('Units','normalized','Position',[0.13 0.16 0.55 0.82]);

handle0 = plot(datTable.Tc(cond2),datTable.ln_CH4_H2O(cond2),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 1 1]-0.2,'MarkerSize',7,...
    'LineWidth',1);
hold on
handle1 = plot(Tc,yM06L,'-k','LineWidth',2);
hold on
handle2 = plot(Tc,yHCTH,'--k','LineWidth',2);
hold on
handle3 = plot(Tc,yHC95_1,'LineWidth',2,'Color',[0.984 0.604 0.600]);
hold on
handle4 = plot(Tc,yHC95_2,'LineWidth',2,'Color',[0.890 0.102 0.110]);
hold on
handle5 = plot(Tc,Horibe95Calcklna-Cerrai54Calcklna-Horita94_h2o_l_g_klna,'LineWidth',2,'Color',[1.000 0.498 0.000]);
hold on
handle6 = plot(Tc,-Richet77cklna-Horita94_h2o_l_g_klna,'LineWidth',2,'Color',[0.122 0.471 0.706]);

ylim([-330 -120])
xlim([-10 60])

xlabel('T (\circC)')
ylabel(['1000ln^2\alpha_{CH_4-H_2O(l)} (' char(8240) ')'])
legend([handle0 handle3 handle4 handle5 handle6 handle1 handle2],...
    {'Environmental','HC95+R76','HC95+S49',...
    'HC95+C54','R77','M06-L','HCTH'},'Location','southeast')
set(gca,'FontSize',14,'YTick',linspace(-320,-120,6))
title('a','Units','normalized','Position',[0.05 0.92 0],'FontSize',15)

axes('Units','normalized','Position',[0.755 0.65 0.22 0.33])
n = 65;
h1 = histogram(datTable.ln_CH4_H2O(cond2)-yHCTH_data',n,...
    'FaceColor','r','LineStyle','none');
hold on
h2 = histogram(datTable.ln_CH4_H2O(cond2)-yM06L_data',n,...
    'FaceColor','k','LineStyle','none');
hold on
line([0 0],[0 40],'LineWidth',2,'Color','k')
xlim([-50 150])
legend([h2,h1],{'M06-L','HCTH'},'Location','northeast','box','off')
% xlabel('\Delta1000ln^2\alpha(obs-exp)')
xlabel('\Delta1000ln^2\alpha_{CH_4-H_2O}')
ylabel('n')
set(gca,'FontSize',12)
title('b','Units','normalized','Position',[0.1 0.8 0],'FontSize',15)

axes('Units','normalized','Position',[0.755 0.16 0.22 0.33])
n = 30;
histogram(datTable.ln_CO2_CH4(cond2)-CHdat_kln13aeq_theo(cond2)',n,...
    'FaceColor','k','LineStyle','none');
hold on
line([0 0],[0 30],'LineWidth',2,'Color','k')
% xlim([-50 100])
% xlabel('\Delta1000ln^{13}\alpha(obs-exp)')
xlabel('\Delta1000ln^{13}\alpha_{CO_2-CH_4}')
ylabel('n')
set(gca,'FontSize',12)
title('c','Units','normalized','Position',[0.1 0.8 0],'FontSize',15)

local_path = '/Users/jonathag/Dropbox (Weizmann Institute)/Apps/Overleaf/EFF paper/figures/';
print('-f2',[local_path,'Fig_7'],'-dpng','-r400');
print('-f2',[local_path,'Fig_7'],'-depsc');

%%
% figure
% axes('Units','normalized','Position',[0.1 0.1 0.7 0.8])
% plot(datTable.Tc(cond2),datTable.ln_CH4_H2O(cond2)-yM06L_data','o',...
%     'MarkerEdgeColor','k','MarkerFaceColor','w')
% hold on
% line([-5 60],[0 0],'LineWidth',2,'Color','k')
% ylim([-50 50])
% xlim([-5 55])
% set(gca,'FontSize',14)
% 
% axes('Units','normalized','Position',[0.81 0.1 0.15 0.8])
% histogram(datTable.ln_CH4_H2O(cond2)-yM06L_data',90,...
%     'FaceColor','k','LineStyle','none');
% set(gca,'XDir','reverse','XTick',[])
% camroll(270)
% xlim([-50 50])
% set(gca,'FontSize',14)
% 

