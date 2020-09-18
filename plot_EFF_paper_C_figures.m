% 5.2.20
% Methanogenesis EFF paper
% Plot beta values for carbon and hydrogen isotopes

% Figures to plot
fig(1) = 1;

file_name    = 'alpha-beta.xlsx';
betaCdat     = readtable(file_name,'Sheet',1);
gammaCHdat   = readtable(file_name,'Sheet',6);
gammaHHdat   = readtable(file_name,'Sheet',8);
clumpedCHdat = readtable(file_name,'Sheet',5);
clumpedHHdat = readtable(file_name,'Sheet',7);

plot_EFF_Cfig(betaCdat,gammaCHdat,clumpedCHdat,...
    gammaHHdat,clumpedHHdat,fig,sav)

if close_val == 1; close all; end

function [] = ...
    plot_EFF_Cfig(betCdat,gammaCHdat,clumpedCHdat,...
    gammaHHdat,clumpedHHdat,fig,sav)
pathstr = '/Users/jonathag/Dropbox (Weizmann Institute)/Jonathan_with_Itay/Paper 2 EFFs/Figures/';
addpath(genpath('/Users/jonathag/Dropbox (Weizmann Institute)/Jonathan_with_Itay/Code/legendflex-pkg-master'))

figWidth     = 14;
figHeight    = 18;
figFontSize  = 15;
figLineWidth = 3;

% Tc = betC(:,1);
% Tk = Tc + 273.15;
Tk = betCdat.temp_C + 273.15;
% gas phase
betaC(1,:) = betCdat.co2_g; % CO2 (g)
betaC(2,:) = betCdat.ch4_g; % CH4 (g)

% aqueous phase
betaC(3,:)  = betCdat.co2;      % CO2
betaC(4,:)  = betCdat.chomfr;   % CHO-MFR
betaC(5,:)  = betCdat.choh4mpt; % CHO-H4MPT
betaC(6,:)  = betCdat.chh4mpt;  % CH-H4MPT
betaC(7,:)  = betCdat.ch2h4mpt; % CH2-H4MPT
betaC(8,:)  = betCdat.ch3h4mpt; % CH3-H4MPT
betaC(9,:)  = betCdat.ch3scom;  % CH3-S-CoM
betaC(10,:) = betCdat.ch4;      % CH4

betaC(11,:) = betCdat.methanol;      % Methanol
betaC(12,:) = betCdat.acetate_CH3;   % Acetate (-CH3)
betaC(13,:) = betCdat.acetylcoA_CH3; % Acetyl-CoA (-CH3)

for i = 1:5
    alphaC(i+1,:) = betaC(i+3,:)./betaC(i+4,:);
end
alphaC(1,:)  = betaC(1,:)./betaC(4,:);   % CO2(g) --> CHO-MFR
alphaC(7,:)  = betaC(9,:)./betaC(2,:);   % CH3-S-CoM --> CH4(g)
alphaC(8,:)  = betaC(11,:)./betaC(9,:);  % Methanol --> CH3-S-CoM
alphaC(9,:)  = betaC(12,:)./betaC(13,:); % Acetate --> Acetyl-CoA
alphaC(10,:) = betaC(13,:)./betaC(8,:);  % Acetyl-CoA --> CH3-H4MPT
% Net EFFs
alphaC(11,:) = betaC(1,:)./betaC(2,:);   % CO2(g) --> CH4(g)
alphaC(12,:) = betaC(12,:)./betaC(2,:);  % Acetate --> CH4(g)
alphaC(13,:) = betaC(11,:)./betaC(2,:);  % Methanol --> CH4(g)

% EFFs between CO2 and CHO-MFR, CHO-H4MPT, CH-H4MPT, CH2-H4MPT, CH3-H4MPT,
% CH3-S-CoM, CH4, Methanol, Acetate (-CH3), Acetyl-CoA (-CH3)
alpha_CO2  = betaC(1,:)./betaC(4:13,:);

betaLoc = betaC(:,1);
%% Expeimental data from literature

% Horita 2001 - CO2-CH4
Horita01Tc   = [200.6 251.2 301.4 348.5 409 449 496.5 550 598.5];
Horita01klna = [34.24 29.33 25.19 21.98 19.26 16.93 15.29 13.87 12.7];

% Kueter 2019 - CO2-CH4
Kueter19Tc   = [300 400 600 800 1000 1200];
Kueter19klna = [25.2 20.1 11.7 8.3 6.9 4.7];

%% Theoretical data from the literture
Bottinga69klna = 2.28e6./(Tk.^2) + 15.176e3./Tk - 8.38; % CO2-CH4, 0 - 700 oc
Richet77klna   = -0.62e9./(Tk.^3) + 6.616e6./(Tk.^2) + 6.04e3./Tk - 3.08; % CO2-CH4, 400 - 1300 oc

%%
col(1,:)  = [0.655 0.808 0.890]; % Fmd / CO2
col(2,:)  = [0.122 0.471 0.706]; % Ftr / CHO-MFR
col(3,:)  = [0.698 0.875 0.541]; % Mch / CHO-H4MPT
col(4,:)  = [0.200 0.627 0.173]; % Mtd / CH-H4MPT
col(5,:)  = [0.984 0.604 0.600]; % Mer / CH2-H4MPT
col(6,:)  = [0.890 0.102 0.110]; % Mtr / CH3-H4MPT
col(7,:)  = [0.992 0.749 0.435]; % Mcr / CH3-S-CoM
col(8,:)  = [1.000 0.498 0.000]; % Mta / CH4
col(9,:)  = [0.792 0.698 0.839]; % Acs / CH3-OH
col(10,:) = [0.416 0.239 0.604]; % Cdh / CH3-COOH
col(11,:) = [0.694 0.349 0.157]; % 

%% 13BETA VALUES FOR 1-100 oc
if fig(1) == 1
    textLab1 = {'CO_2(g) [+4]','CHO-MFR [+2]','CHO-H_4MPT [+2]','CH-H_4MPT [+2]','CH_2-H_4MPT [0]',...
        'CH_3-H_4MPT [-2]','CH_3-SCoM [-2]','CH_3-OH [-2]','CH_3-COOH [-3]','CH_4(g) [-4]'};
    
    x1 = 1e3./Tk;
    x1a = 1e3./(0+273.15);
    
    figure(1)
    set(1,'Units','centimeters','Position',[10 5 figWidth figHeight])
    pos = [1,4:9,11:12,2];
    pos_col = [1:7 9:10 8];
    interp = 'tex';
    for ja = 1:length(pos)
        plot(x1(1:51),betaC((pos(ja)),1:51),'Color',col(pos_col(ja),:),...
            'LineWidth',figLineWidth)
        if pos(ja) == 4
            t = text(x1a+0.02,betaLoc(pos(ja),1)+0.005,textLab1{ja},'Interpreter',interp);
        elseif pos(ja) == 5
            t = text(x1a+0.02,betaLoc(pos(ja),1)-0.002,textLab1{ja},'Interpreter',interp);
        elseif pos(ja) == 11
            t = text(x1a+0.02,betaLoc(pos(ja),1)+0.006,textLab1{ja},'Interpreter',interp);
        elseif pos(ja) == 12
            t = text(x1a+0.02,betaLoc(pos(ja),1)-0.002,textLab1{ja},'Interpreter',interp);
        elseif pos(ja) == 2
            t = text(x1a+0.02,betaLoc(pos(ja),1)-0.002,textLab1{ja},'Interpreter',interp);
        else
            t = text(x1a+0.02,betaLoc(pos(ja),1),textLab1{ja},'Interpreter',interp);
        end
        t.FontName = 'Helvetica';
        t.FontSize = 14;
        hold on
    end
    box off
    xlabel('1000/T (K)')
    ylabel('^{13}\beta')
    xlim([2.6 4.35])
    ylim([1.08 1.23])
    grid on
    grid minor
    set(gca,'FontSize',figFontSize)
    if sav == 1
        print('-f1',sprintf('%s13beta',pathstr),'-dpng','-r600');
        print('-f1',sprintf('%s13beta',pathstr),'-depsc');
    end
end
%% 13ALPHA VALUES FOR ENZYMES
if fig(2) == 1
    textLab2 = {'Fmd';'Ftr';'Mch';'Mtd';'Mer';'Mtr';'Mcr';'Mta';'Acs';'Cdh'};
    alphaLoc = alphaC(:,51);
    x1 = betCdat.temp_C;
    
    figure(2)
    set(2,'Units','centimeters','Position',[10 5 figWidth figHeight])
    for ja = 1:10
        ax1(ja) = plot(x1(1:51),1000.*log(alphaC((ja),1:51)),'Color',col(ja,:),...
            'LineWidth',figLineWidth);
            t = text(103,1000.*log(alphaLoc(ja,1)),textLab2{ja});
        t.FontName = 'Helvetica';
        t.FontSize = 16;
        hold on
    end
    box off
    xlabel('T (\circC)')
    ylabel(['1000ln^{13}\alpha^{eq} (' char(8240) ')'])
    xlim([0 120])
    ylim([-5 22])
    grid on
    grid minor
    set(gca,'FontSize',figFontSize)
    if sav == 1
        print('-f2',sprintf('%s13alpha',pathstr),'-dpng','-r600');
        print('-f2',sprintf('%s13alpha',pathstr),'-depsc');
    end
end

%% 13ALPHA VALUES RELATIVE TO CO2
if fig(3) == 1
    textLab2 = {'CHO-MFR [+2]','CHO-H_4MPT [+2]','CH-H_4MPT [+2]',...
        'CH_2-H_4MPT [0]','CH_3-H_4MPT [-2]','CH_3-SCoM [-2]',...
        'CH_4(g) [-4]','CH_3-OH [-2]','CH_3-COOH [-3]','CH_3-CoA [-3]'};
    alphaLoc = alpha_CO2(:,51);
    x1 = betCdat.temp_C;
    
    figure(3)
    set(3,'Units','centimeters','Position',[10 5 figWidth figHeight])
    for ja = 1:9
        ax1(ja) = plot(x1(1:51),1000.*log(alpha_CO2((ja),1:51)),...
            'Color',col(ja+1,:),...
            'LineWidth',figLineWidth);
        if ja == 3
            t = text(103,1000.*log(alphaLoc(ja,1))-0.5,textLab2{ja});
        elseif ja == 8
            t = text(103,1000.*log(alphaLoc(ja,1))-1,textLab2{ja});
        else
            t = text(103,1000.*log(alphaLoc(ja,1)),textLab2{ja});
        end
        t.FontName = 'Helvetica';
        t.FontSize = 13;
        hold on
    end
    xlabel('T (\circC)')
    ylabel(['1000ln^{13}\alpha^{eq}_{CO_2' char(8211) 'X} (' char(8240) ')'])
    xlim([0 160])
    grid on
    grid minor
    box off
    set(gca,'FontSize',figFontSize)
    if sav == 1
        print('-f3',sprintf('%s13alpha_co2',pathstr),'-dpng','-r600');
        print('-f3',sprintf('%s13alpha_co2',pathstr),'-depsc');
    end
end

%% COMPARISON OF 13ALPHA VALUES FOR CO2-CH4
if fig(4) == 1
    figure(4)
    set(4,'Units','centimeters','Position',[10 5 figWidth figHeight])
    x_axis = Tk-273.15;
    x_axis = 1000./Tk;
    
    h3 = plot(x_axis,1000.*log(betaC(1,:)./betaC(2,:)),'k','LineWidth',figLineWidth);
    hold on

    h2 = plot(x_axis,Richet77klna,'k','LineWidth',figLineWidth);
    h2.Color = [31 120 180]./255;
    h2.LineStyle = ':';
    hold on
    
    h4 = plot(1000./(Horita01Tc+273.15),Horita01klna,'o');
    h4.MarkerEdgeColor = 'k';
    h4.MarkerFaceColor = 'w';
    h4.MarkerSize = 10;
    h4.LineWidth = 0.5;
    hold on
    
    h5 = plot(1000./(Kueter19Tc(1:3)+273.15),Kueter19klna(1:3),'o');
    h5.MarkerEdgeColor = 'k';
    h5.MarkerFaceColor = [1 1 1]-0.5;
    h5.MarkerSize = 10;
    h5.LineWidth = 0.5;
    hold on
    
    l = legend([h4 h5 h2 h3],...
        {'Horita (2001)',...
        'Kueter et al. (2019)',...
        'Richet et al. (1977)',...
        'This work'});
    l.Location = 'northwest';
    l.Box = 'on';
    l.FontSize = 14;
    grid on
    xlabel('1000/T (K)')
    ylabel(['1000ln^{13}\alpha_{CO_2(g)-CH_4(g)} (' char(8240) ')'])
    set(gca,'FontSize',figFontSize')
    xlim([1 3.8])
    ylim([5 80])
    if sav == 1
        print('-f4',sprintf('%s13alpha_co2_ch4_comp',pathstr),'-dpng','-r600');
        print('-f4',sprintf('%s13alpha_co2_ch4_comp',pathstr),'-depsc');
    end
end

%% CLUMPED ISOTOPOLOGUES (13--D)
if fig(5) == 1
    figure(5)
    set(5,'Units','centimeters','Position',[10 5 figWidth figHeight])
%     Deltalab = {'CHO-MFR','CHO-H$_4$MPT','CH-H$_4$MPT','CH$_2$-H$_4$MPT',...
%         'CH$_3$-H$_4$MPT','CH$_3$-S-CoM','CH$_4$','CH$_3$-OH','CH$_3$-COOH','CH$_3$-CSCoA'}; 
    Deltalab = {'13CDO-MFR','13CDO-H4MPT','13CD-H4MPT','13CHD-H4MPT',...
        '13CH2D-H4MPT','13CH2D-SCoM','13CH3D','13CH2D-OH',...
        '13CH2D-COOH','13CH2D-CSCoA'}; 
    for i = 1:10
        ax(i) = plot(clumpedCHdat.T_C(1:2:51),...
                     clumpedCHdat{1:2:51,i+1},'Color',col(i+1,:),...
            'LineWidth',figLineWidth);
        hold on
        t = text(103,clumpedCHdat{51,i+1},Deltalab{i});
        t.FontName = 'Helvetica';
        t.FontSize = 14;
    end
    xlim([0 150])
    ylim([2.7 6.7])
    grid on
    grid minor
    xlabel('T (\circC)')
    ylabel(['\Delta_i^{eq} (' char(8240) ')'])
    set(gca,'FontSize',figFontSize)
    box off
    if sav == 1
        print('-f5',sprintf('%sDelta_132',pathstr),'-dpng','-r600');
        print('-f5',sprintf('%sDelta_132',pathstr),'-depsc');
    end
end

%% CLUMPED ISOTOPOLOGUES (D--D)
if fig(7) == 1
    figure(7)
    set(7,'Units','centimeters','Position',[10 5 figWidth figHeight])
    Delta_DD_lab = {'12CD2-H4MPT','12CHD2-H4MPT','12CHD2-SCoM',...
        '12CH2D2','12CHD2-COO','12CHD2-COSCoA','12CHD2-OH'}; 
    
    for i = 1:7
        ax(i) = plot(clumpedHHdat.T_C(1:2:51),...
                     clumpedHHdat{1:2:51,i+1},'Color',col(i+4,:),...
            'LineWidth',figLineWidth);
        hold on
        t = text(103,clumpedHHdat{51,i+1},Delta_DD_lab{i});
        t.FontName = 'Helvetica';
        t.FontSize = 14;
    end
    xlim([0 150])
    ylim([7.5 23])
    grid on
    grid minor
    xlabel('T (\circC)')
    ylabel(['\Delta_i^{eq} (' char(8240) ')'])
    set(gca,'FontSize',figFontSize)
    box off
    if sav == 1
        print('-f7',sprintf('%sDelta_22',pathstr),'-dpng','-r600');
        print('-f7',sprintf('%sDelta_22',pathstr),'-depsc');
    end
end

%% COMPARISON OF DELTA(12CH2D2) OF METHANE
if fig(8) == 1
    figure(8)
    set(8,'Units','centimeters','Position',[10 5 figWidth figHeight])
    
    x = Tk - 273.15;
    
    eldridge2019_12CH2D2 = -9.67634e15./Tk.^6 + 171.917e12./Tk.^5 - ...
        1248.19e9./Tk.^4 + 4302.83e6./Tk.^3 - 4486.6e3./Tk.^2 + 1862.58./Tk;
    young2017_12CH2D2 = 1000.*log(1 + 0.183798./Tk - 785.483./Tk.^2 + ...
        1056.280e3./Tk.^3 + 93.7307e6./Tk.^4 - ...
        89.1948e9./Tk.^5 + 9.90173e12./Tk.^6);
    
    eldridge2019_exp_Tc = [1.20 25.00 50.50 75.70 127.80 165.40 250.00 300.00 350.00 400.00 500.00];
    eldridge2019_exp_12CH2D2 = [23.50 19.64 14.09 14.05 9.64 8.88 5.00 1.86 3.44 3.46 -0.06];
    eldridge2019_exp_SE = [1.10 0.67 0.33 1.40 0.42 1.02 1.44 1.22 1.21 0.28 0.64];
        
    han2 = plot(x,young2017_12CH2D2,'Color',[31 120 180]./255,'LineWidth',figLineWidth);
    hold on
    han1_g = plot(clumpedHHdat.T_C,clumpedHHdat.CH2D2,'Color','k',...
            'LineWidth',figLineWidth);
    hold on
    han1_SMD = plot(clumpedHHdat.T_C,clumpedHHdat.CH2D2_SMD,'Color','k',...
            'LineWidth',figLineWidth,'LineStyle','--');
    hold on
    han3 = plot(x,eldridge2019_12CH2D2,'Color',[227 26 28]./255,...
        'LineStyle','--','LineWidth',figLineWidth);
    hold on
    han4 = errorbar(eldridge2019_exp_Tc,eldridge2019_exp_12CH2D2,eldridge2019_exp_SE,...
        'LineStyle','none','Marker','o','CapSize',4,'Color',[227 26 28]./255,...
        'MarkerFaceColor','w','MarkerSize',9,'LineWidth',1.5);
    grid on
    grid minor
    legend([han2 han3 han4 han1_g han1_SMD],...
        {'Young et al. (2017)','Eldridge et al. (2019)',...
        'Eldridge et al. (2019)','This work (M06L)','This work (M06L-SMD)'},...
        'Location','northeast','Box','off');
    xlabel('T (\circC)')
    ylabel(['\Delta_{^{12}CH_2D_2}^{eq} (' char(8240) ')'])
    ylim([-1.5 25])
    xlim([0 700])
    set(gca,'FontSize',figFontSize)
    box off
    if sav == 1
%         print('-f8',sprintf('%sDelta_12CH2D2',pathstr),'-dpng','-r600');
        print('-f8',sprintf('%sDelta_12CH2D2',pathstr),'-depsc');
    end
end

%% COMPARISON OF DELTA(12CH3D) OF METHANE
if fig(9) == 1
    figure(9)
    set(9,'Units','centimeters','Position',[10 5 figWidth figHeight])
    
    x = Tk - 273.15;
    
    eldridge2019_13CH3D = 1.47348e19./Tk.^7 - 2.08648e17./Tk.^6 + ...
        1.19810e15./Tk.^5 - 3.54757e12./Tk.^4 + 5.54476e9./Tk.^3 - ...
        3.49294e6./Tk.^2 + 8.89370e2./Tk;
    
    young2017_13CH3D = 1000.*log(1 + 0.0355502./Tk - 433.038./Tk.^2 + ...
        1270210./Tk.^3 - 5.94804e8./Tk.^4 + 1.196630e11./Tk.^5 - 9.07230e12./Tk.^6);
    
    webb_miller_2014_Tk = [300 400 500 600];
    webb_miller_2014_13CH3D_PIMC = [5.73 3.5 2.29 1.51];
    
    eldridge2019_exp_Tc = [1.20 25.00 50.50 75.70 127.80 165.40 250.00 300.00 350.00 400.00 500.00];
    eldridge2019_exp_13CH3D_i = [7.00 5.62 5.07 4.47 3.35 2.93 2.27 1.26 1.46 1.44 0.76];
    eldridge2019_exp_13CH3D = (exp(eldridge2019_exp_13CH3D_i./1000)-1).*1000;
    eldridge2019_exp_SE = [0.06 0.05 0.17 0.15 0.08 0.21 0.20 0.11 0.10 0.36 0.06];
    
    han1 = plot(clumpedCHdat.T_C,clumpedCHdat.ch4,'Color','k',...
            'LineWidth',figLineWidth);
    hold on
    han2 = plot(x,young2017_13CH3D,'Color',[31 120 180]./255,...
        'LineWidth',figLineWidth,'LineStyle',':');
    hold on
    han3 = plot(x,eldridge2019_13CH3D,'Color',[227 26 28]./255,...
        'LineStyle','--','LineWidth',figLineWidth);
    hold on
    han4 = plot(webb_miller_2014_Tk-273.15,webb_miller_2014_13CH3D_PIMC,...
        'o','MarkerFaceColor','w','MarkerEdgeColor',[51 160 44]./255,...
        'MarkerSize',9,'LineWidth',1.5);
    hold on
    han5 = errorbar(eldridge2019_exp_Tc,eldridge2019_exp_13CH3D,eldridge2019_exp_SE,...
        'LineStyle','none','Marker','o','CapSize',4,'Color',[227 26 28]./255,...
        'MarkerFaceColor','w','MarkerSize',9,'LineWidth',1.5);
    grid on
    grid minor
    legend([han2 han3 han4 han5 han1],...
        {'Young et al. (2017)','Eldridge et al. (2019)',...
        'Webb & Miller (2014)','Eldridge et al. (2019)','This work'},...
        'Location','northeast','Box','off');
    xlabel('T (\circC)')
    ylabel(['\Delta_{^{13}CH_3D}^{eq} (' char(8240) ')'])
    xlim([0 700])
    set(gca,'FontSize',figFontSize)
    box off
    if sav == 1
        print('-f9',sprintf('%sDelta_13CH3D',pathstr),'-depsc');
    end
end

end