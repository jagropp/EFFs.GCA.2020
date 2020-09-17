% 19.11.18
% EFF paper
% Plot beta values for carbon and hydrogen isotopes
clear

local_path = '/Users/jonathag/Dropbox (Weizmann Institute)/Paper 2 Methanogenesis EFFs/';
file_name = 'alpha-beta.xlsx';
betaCdat  = readtable([local_path,file_name],'Sheet',1);
betaCHdat = readtable([local_path,file_name],'Sheet',5);
gammaCHdat = readtable([local_path,file_name],'Sheet',6);
alphaCHdat = readtable([local_path,file_name],'Sheet',4);

fig = [1 1 1 1 1];
sav  = 0;
plot_EFF_Cfig(betaCdat,betaCHdat,alphaCHdat,gammaCHdat,fig,sav)
% close all

function [] = plot_EFF_Cfig(betCdat,betCHdat,alphaCHdat,gammaCHdat,fig,sav)
pathstr = '/Users/jonathag/Dropbox (Weizmann Institute)/Paper 2 Methanogenesis EFFs/Figures/';
addpath(genpath('/Users/jonathag/Dropbox (Weizmann Institute)/Jonathan_(with_Itay)/Code/legendflex-pkg-master'))

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

betaLoc = betaC(:,1);
%% Expeimental data from literature

% Horita 2001 - CO2-CH4
Horita01Tc   = [200.6 251.2 301.4 348.5 409 449 496.5 550 598.5];
Horita01klna = [34.24 29.33 25.19 21.98 19.26 16.93 15.29 13.87 12.7];
Horita01std  = [0.12 0.19 0.11 0.21 0.11 0.1 0.14 0.19 0.15];

% Kueter 2019 - CO2-CH4
Kueter19Tc   = [300 400 600 800 1000 1200];
Kueter19klna = [25.2 20.1 11.7 8.3 6.9 4.7];
Kueter19std  = [0.5 0.6 0.7 0.2 0.2 0.2];

%% Theoretical data from the literture
Bottinga69klna = 2.28e6./(Tk.^2) + 15.176e3./Tk - 8.38; % CO2-CH4, 0 - 700 oc
Richet77klna   = -0.62e9./(Tk.^3) + 6.616e6./(Tk.^2) + 6.04e3./Tk - 3.08; % CO2-CH4, 400 - 1300 oc
Mook74klna1    = 9.866e3./Tk - 24.12; % [HCO3-] - CO2(aq) 25-125 oc
Mook74klna2    = 9.552e3./Tk - 24.10; % [HCO3-] - CO2(aq) 25-125 oc
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

%% FIGURE 1 - 13BETA VALUES FOR 1-100 oc
if fig(1) == 1
    textLab1 = {'CO_2(g) [+4]','CH_4(g) [-4]','CHO-MFR [+2]',...
        'CHO-H_4MPT [+2]','CH-H_4MPT [+2]','CH_2-H_4MPT [0]',...
        'CH_3-H_4MPT [-2]','CH_3-S-CoM [-2]','CH_3-OH [-2]','CH_3-COOH [-3]'};
    textLab1 = {'CO2(g) [+4]','CHO-MFR [+2]','CHO-H4MPT [+2]','CH-H4MPT [+2]','CH2-H4MPT [0]',...
        'CH3-H4MPT [-2]','CH3-S-CoM [-2]','CH3-OH [-2]','CH3-COOH [-3]','CH4(g) [-4]'};
    
    x1 = 1e3./Tk;
    x1a = 1e3./(0+273.15);
    
    figure(1)
    set(1,'Units','centimeters','Position',[10 5 9 9])
    pos = [1,4:9,11:12,2];
    for ja = 1:length(pos)
        plot(x1(1:51),betaC((pos(ja)),1:51),'Color',col(ja,:))
        if pos(ja) == 4
            t = text(x1a+0.02,betaLoc(pos(ja),1)+0.005,textLab1{ja});
        elseif pos(ja) == 5
            t = text(x1a+0.02,betaLoc(pos(ja),1)-0.002,textLab1{ja});
        elseif pos(ja) == 11
            t = text(x1a+0.02,betaLoc(pos(ja),1)+0.006,textLab1{ja});
        elseif pos(ja) == 12
            t = text(x1a+0.02,betaLoc(pos(ja),1)-0.002,textLab1{ja});
        elseif pos(ja) == 2
            t = text(x1a+0.02,betaLoc(pos(ja),1)-0.002,textLab1{ja});
        else
            t = text(x1a+0.02,betaLoc(pos(ja),1),textLab1{ja});
        end
        t.FontName = 'Helvetica';
        t.FontSize = 8;
        hold on
    end
    box off
    xlabel('1000/T (K)')
    ylabel(['^{13}\beta (' char(8240) ')'])
    xlim([2.6 4.45])
    % ylim([0.98 1.09])
    grid on
    grid minor
    if sav == 1
        print('-f1',sprintf('%s13beta',pathstr),'-dpng','-r600');
        print('-f1',sprintf('%s13beta',pathstr),'-depsc');
    end
end
%% FIGURE 2 - 13ALPHA VALUES
if fig(2) == 1
    textLab2 = {'Fmd';'Ftr';'Mch';'Mtd';'Mer';'Mtr';'Mcr';'Mta';'Acs';'Cdh'};
    alphaLoc = alphaC(:,51);
    x1 = betCdat.temp_C;
    
    offset = 0.4;
    figure(2)
    set(2,'Units','centimeters','Position',[10 5 9 9])
    for ja = 1:10
        ax1(ja) = plot(x1(1:51),1000.*log(alphaC((ja),1:51)),'Color',col(ja,:));
%         if ja == 4
%             t = text(103,1000.*log(alphaLoc(ja,1))+offset,textLab2{ja});
%         elseif ja == 6
%             t = text(103,1000.*log(alphaLoc(ja,1))-offset,textLab2{ja});
%         elseif ja == 10
%             t = text(103,1000.*log(alphaLoc(ja,1))+offset,textLab2{ja});
%         elseif ja == 3
%             t = text(103,1000.*log(alphaLoc(ja,1))-offset,textLab2{ja});
%         elseif ja == 1
%             t = text(103,1000.*log(alphaLoc(ja,1))-0.4,textLab2{ja});
%         elseif ja == 9
%             t = text(103,1000.*log(alphaLoc(ja,1))+0.4,textLab2{ja});
%         else
%             t = text(103,1000.*log(alphaLoc(ja,1)),textLab2{ja});
%         end
%         t.FontName = 'Helvetica';
        hold on
    end
    hl(1).leg = legendflex(ax1, textLab2, 'anchor', [3 3], ...
        'buffer', [-3 -3], ...
        'ncol', 2, ...
        'fontsize', 8, ...
        'xscale', 1, ...
        'box', 'on',...
        'FontName','Helvetica');
    box off
    xlabel('T (\circC)')
    ylabel(['1000ln^{13}\alpha (' char(8240) ')'])
    xlim([0 100])
    ylim([-5 28])
    grid on
    grid minor
    if sav == 1
        print('-f2',sprintf('%s13alpha',pathstr),'-dpng','-r600');
        print('-f2',sprintf('%s13alpha',pathstr),'-depsc');
    end
end
%% FIGURE 3 - CO2-CH4 COMPARISON
if fig(3) == 1
    figure(3)
    set(3,'Units','centimeters','Position',[10 5 9 9])
    h3 = plot(1000./Tk,1000.*log(betaC(1,:)./betaC(2,:)),'k');
    hold on
    h1 = plot(1000./Tk,Bottinga69klna);
    h1.Color = [0.651 0.808 0.890];%[1 1 1].*0.6;
    h1.LineStyle = ':';
    hold on
    h2 = plot(1000./Tk,Richet77klna,'k');
    h2.Color = [0.122 0.471 0.706];%[1 1 1].*0.4;
    h2.LineStyle = ':';
    hold on
    h4 = plot(1000./(Horita01Tc+273.15),Horita01klna,'o');
    h4.MarkerEdgeColor = 'k';
    h4.MarkerFaceColor = [0.200 0.627 0.173];%[1 1 1].*0.9;
    h4.LineWidth = 0.5;
    hold on
    h5 = plot(1000./(Kueter19Tc(1:3)+273.15),Kueter19klna(1:3),'o');
    h5.MarkerEdgeColor = 'k';
    h5.MarkerFaceColor = [0.698 0.875 0.541];%[1 1 1].*0.6;
    h5.LineWidth = 0.5;
    hold on
    l = legend([h4 h5 h1 h2 h3],...
        {'Horita (2001)',...
        'Kueter et al. (2019)',...
        'Bottinga (1969)',...
        'Richet et al. (1977)',...
        'This work'});
    l.Location = 'southeast';
    l.Box = 'on';
    % l.EdgeColor = 'w';
    grid on
    grid minor
    xlabel('1000/T (K)')
    ylabel(['1000ln^{13}\alpha_{CO_2(g)-CH_4(g)} (' char(8240) ')'])
    if sav == 1
        print('-f3',sprintf('%s13alpha_co2_ch4_comp',pathstr),'-dpng','-r600');
        print('-f3',sprintf('%s13alpha_co2_ch4_comp',pathstr),'-depsc');
    end
end

%% FIGURE 4 - Clumped isotopologues
if fig(4) == 1
    figure(4)
    set(4,'Units','centimeters','Position',[10 5 9 9])
    Deltalab = {'CHO-MFR','CHO-H_4MPT','CH-H_4MPT','CH_2-H_4MPT','CH_3-H_4MPT',...
         'CH_3-S-CoM','CH_4','CH_3-OH','CH_3-COOH','CH_3-CSCoA'}; 
     Deltalab = {'CHO-MFR','CHO-H4MPT','CH-H4MPT','CH2-H4MPT','CH3-H4MPT',...
         'CH3-S-CoM','CH4','CH3-OH','CH3-COOH','CH3-CSCoA'}; 
    for i = 1:10
        ax(i) = plot(betCHdat.T_C(1:51),betCHdat{1:51,i+1},'Color',col(i,:));
        hold on
%         t = text(103,betCHdat{51,i+1},Deltalab{i});
%         t.FontName = 'Helvetica';
    end
    hl(1).leg = legendflex(ax, Deltalab, 'anchor', [3 3], ...
        'buffer', [-3 -3], ...
        'ncol', 2, ...
        'fontsize', 8, ...
        'xscale', 0.7, ...
        'box', 'on',...
        'FontName','Helvetica');
    xlim([0 100])
    ylim([2.6 7.5])
    grid on
    grid minor
    xlabel('T (\circC)')
    ylabel(['\Delta_i (' char(8240) ')'])
    if sav == 1
        print('-f4',sprintf('%sDelta',pathstr),'-dpng','-r600');
        print('-f4',sprintf('%sDelta',pathstr),'-depsc');
    end
end

%% FIGURE 5 - Clumped isotopologues
if fig(5) == 1    
    figure(5)
    set(5,'Units','centimeters','Position',[10 5 9 9])
    alphalab_leg = {'Fmd','Ftr','Mch','Mtd','Mer','Mtr','Mcr','Mta','Acs','Cdh'}; 
  
    pos_s = [2:5,7:8,13:15];
    pos_p = [1,9,11];
    pos_s_col = [2:6,8:10];
    pos_p_col = [1,4,5];
    
    maxval = 51;
    x = gammaCHdat.T_C(1:maxval);
    hold all
    ax2(1)  = plot(x,10.*gammaCHdat.Fmd_p(1:maxval),'Color',col(1,:));
    ax2(2)  = plot(x,gammaCHdat.Ftr_s(1:maxval),    'Color',col(2,:));
    ax2(3)  = plot(x,gammaCHdat.Mch_s(1:maxval),    'Color',col(3,:));
    ax2(4)  = plot(x,gammaCHdat.Mtd_s(1:maxval),    'Color',col(4,:));
    ax2(5)  = plot(x,gammaCHdat.Mer_R_s(1:maxval),  'Color',col(5,:));
    ax2(6)  = plot(x,gammaCHdat.Mtr_s(1:maxval),    'Color',col(6,:));
    ax2(7)  = plot(x,gammaCHdat.Mcr_s(1:maxval),    'Color',col(7,:));
    ax2(8)  = plot(x,gammaCHdat.Mta_s(1:maxval),    'Color',col(8,:));
    ax2(9)  = plot(x,gammaCHdat.Acs_s(1:maxval),    'Color',col(9,:));
    ax2(10) = plot(x,gammaCHdat.Cdh_s(1:maxval),    'Color',col(10,:));
    ax2(11) = plot(x,gammaCHdat.Fmd_p(1:maxval),    'Color',col(1,:));
    ax2(12) = plot(x,gammaCHdat.Mtd_p(1:maxval),    'Color',col(4,:));
    ax2(13) = plot(x,gammaCHdat.Mer_p(1:maxval),    'Color',col(5,:));
    ax2(14) = plot(x,gammaCHdat.Mcr_p(1:maxval),    'Color',col(7,:));
    ax2(15) = plot(x,gammaCHdat.Hmd_p(1:maxval),    'Color',col(2,:));
    ax2(16) = plot(x,gammaCHdat.Mer_S_s(1:maxval),  'Color',col(5,:));
    
    for i = 11:15; ax2(i).LineStyle = '-.'; end
    ax2(16).LineStyle = ':'; 

    legendflex(ax2(1:10), alphalab_leg, 'anchor', [5 5], ...
        'buffer', [-10 10], ...
        'ncol', 2, ...
        'fontsize', 8, ...
        'xscale', 0.7, ...
        'box', 'on',...
        'FontName','Helvetica');
    xlim([0 100])
    ylim([0.992 1.0006])
    grid on
    grid minor
    xlabel('1000/T (K)')
    ylabel('\gamma')
    if sav == 1
        print('-f5',sprintf('%s132alpha',pathstr),'-dpng','-r600');
        print('-f5',sprintf('%s132alpha',pathstr),'-depsc');
    end
end
end