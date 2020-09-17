% 19.11.18
% EFF paper
% Plot beta values for hydrogen isotopes

clear
% [betH,betHlab] = xlsread('C:\Jonathan\Dropbox (Weizmann Institute)\Paper 2 Methanogenesis EFFs\alpha-beta.xlsx',4);
% betH(1,:) = [];

local_path = '/Users/jonathag/Dropbox (Weizmann Institute)/Jonathan_with_Itay/Paper 2 EFFs/';
% local_path = '/Users/jonathag/Dropbox (Weizmann Institute)/Apps/Overleaf/EFF paper/';
pathstr    = [local_path,'Figures/']; 
file_name  = 'alpha-beta.xlsx';
betHdat    = readtable([local_path,file_name],'Sheet',2);

fig1 = 0; % 2beta
fig2 = 0; % 1000ln2a
fig3 = 0; % 1000ln2a primary (old)
fig4 = 0; % H2O(l)-H2 comparison
fig5 = 0; % CH4-H2 comparison
fig6 = 0; % CH4-H2O comparison
fig7 = 0; % H2O(l)-H2 comparison with full substituion
fig8 = 0; % CH4-H2 comparison with full substitution
fig9 = 0; % 1000ln2a vs. H2O
sav  = 0;
plot_EFF_Hfig(betHdat,fig1,fig2,fig3,fig4,fig5,...
    fig6,fig7,fig8,fig9,sav,pathstr)
% close all

function [] = plot_EFF_Hfig(betHdat,fig1,fig2,fig3,fig4,...
                            fig5,fig6,fig7,fig8,fig9,sav,pathstr)
marksize = 10;
addpath(genpath('/Users/jonathag/Dropbox (Weizmann Institute)/Jonathan_with_Itay/Code/legendflex-pkg-master'))

figWidth     = 14;
figHeight    = 18;
figFontSize  = 15;
figLineWidth = 3;

Tc = betHdat.T_C;
Tk = betHdat.T_C + 273.15;

% gas phase
betaH(1,:)  = betHdat.HDO_g;      % H2O
betaH(2,:)  = betHdat.C12H3D_g;   % CH4
betaH(3,:)  = betHdat.HD_g;       % H2
% aqueous phase
betaH(4,:)  = betHdat.HDO_l;      % h2o(l)
betaH(5,:)  = betHdat.chomfr;     % formmfr
betaH(6,:)  = betHdat.choh4mpt;   % formh4mpt
betaH(7,:)  = betHdat.chh4mpt;    % menylh4mpt
betaH(8,:)  = betHdat.ch2h4mpt_R; % mlenh4mpt-R
betaH(9,:)  = betHdat.ch2h4mpt_S; % mlenh4mpt-S
betaH(10,:) = betHdat.ch3h4mpt;   % mh4mpt
betaH(11,:) = betHdat.ch3scom;    % mcom
betaH(12,:) = betHdat.f420;       % f420
betaH(13,:) = betHdat.cob;        % cob
betaH(14,:) = betHdat.methanol;   % Methanol
betaH(15,:) = betHdat.aetate;     % Acetate
betaH(16,:) = betHdat.acetylcoa;  % Acetyl-CoA

alphaH_H2O  = betaH([2:3,5:16],:)./betaH(4,:);

% EFFs
alphaH(1,:)  = betaH(4,:)./betaH(5,:);   % Fmd p
alphaH(2,:)  = betaH(5,:)./betaH(6,:);   % Ftr s
alphaH(3,:)  = betaH(6,:)./betaH(7,:);   % Mch s
alphaH(4,:)  = betaH(7,:)./betaH(9,:);   % Mtd s
alphaH(5,:)  = mean(betaH(9:10,:))./betaH(10,:);  % Mer s
alphaH(6,:)  = betaH(10,:)./betaH(11,:); % Mtr s
alphaH(7,:)  = betaH(11,:)./betaH(2,:);  % Mcr s
alphaH(8,:)  = betaH(4,:)./betaH(12,:);  % Frh p
alphaH(9,:)  = betaH(12,:)./betaH(8,:);  % Mtd p
alphaH(10,:) = betaH(12,:)./betaH(10,:); % Mer p
alphaH(11,:) = betaH(4,:)./betaH(13,:);  % Mvh/Hdr p
alphaH(12,:) = betaH(13,:)./betaH(2,:);  % Mcr p
alphaH(13,:) = betaH(14,:)./betaH(11,:); % Mta s
alphaH(14,:) = betaH(15,:)./betaH(16,:); % Acs s
alphaH(15,:) = betaH(16,:)./betaH(10,:); % Cdh s
alphaH(16,:) = betaH(2,:)./betaH(4,:);   % net CH4-H2O(l)
alphaH(17,:) = betaH(2,:)./betaH(3,:);   % net CH4-H2
alphaH(18,:) = betaH(1,:)./betaH(3,:);   % net H2O(g)-H2
alphaH(19,:) = betaH(4,:)./betaH(3,:);   % net H2O(l)-H2
alphaH(20,:) = betaH(14,:)./betaH(2,:);  % net methanol-CH4
alphaH(21,:) = betaH(15,:)./betaH(2,:);  % net acetate-CH4

%% H2O(l)-H2O(g) for 0-374oc
Horita94_h2o_l_g_klna = -0.353e18./Tk.^6 + 42.17e12./Tk.^4 - ...
    309.4e9./Tk.^3 + 963.7e6./Tk.^2 - 1399e3./Tk + 766.2;
Horita94_h2o_l_g_klna(Tc>374) = nan;
%% Expeimental data from literature

% Cerrai 1954 - h2o(g)-h2(g)
Cerrai54Tc = [51.0 59.0 64.0 74.0 77.5 79.4 89.5 97.0 109.0 118.5...
              134.0 148.4 168.0 186.5 202.0 219.0 251.0 280.0 315.0...
              357.0 398.0 445.0 497.0 563.0 609.0 742.0];
Cerrai54klna = [1111.86 1043.80 1018.85 1033.18 1026.04 1000.63 974.56...
                932.16 908.26 900.16 832.91 806.48 765.47 732.37 688.13...
                667.83 615.19 548.12 470.00 438.25 418.71 371.56 329.30...
                307.48 270.03 165.51];
Cerrai54_h2o_l_g = -0.353e18./(Cerrai54Tc+273.15).^6 + 42.17e12./(Cerrai54Tc+273.15).^4 - ...
    309.4e9./(Cerrai54Tc+273.15).^3 + 963.7e6./(Cerrai54Tc+273.15).^2 - 1399e3./(Cerrai54Tc+273.15) + 766.2;
Cerrai54klna_l_g = Cerrai54klna + Cerrai54_h2o_l_g;

% Rolston 1976 - h2o(l)-h2(g)
Rolston76Tc = [6.55 6.55 6.88 12.85 14.85 14.85 18.1 19.24 19.24 19.24...
               22.25 22.25 23.8 23.8 30.87 30.87 32.2 32.2 39.99 40.22...
               40.22 40.84 50.32 51.31 51.31 57.02 60.73 60.73 60.77...
               60.77 60.77 71.3 71.3 81.47 82.13 91.08 91.08 94.76];
Rolston76klna = [1456.99 1447.15 1469.95 1485.01 1415.12 1392.77 1365.07...
                 1386.04 1366.09 1373.97 1359.44 1367.37 1332.10 1306.44...
                 1294.45 1300.74 1294.18 1305.08 1278.71 1259.60 1257.61...
                 1250.47 1203.87 1180.19 1180.19 1176.50 1129.46 1129.46...
                 1153.10 1154.68 1141.03 1087.55 1100.28 1043.80 1040.98...
                 991.77 1019.57 992.14];
             
% Horibe and Craig 1995 - ch4(g)-h2(g)
Horibe1995Tc = [203.3 203.4 217.2 218 220.1 220.1 224 263.3 264.2 264.7...
                300.7 302.2 303.2 304.1 361.8 365.1 365.4 366 403.2 405.5...
                412.3 456.9 457.5 459.9 476.5 499.7 503 505.2 505.7];
Horibe1995klna = [538.049 533.638 509.031 505.718 501.010 501.802 492.893...
                  431.733 431.010 427.835 377.183 377.538 374.710 372.640...
                  303.476 298.386 299.703 297.914 260.636 259.669 256.629...
                  218.680 217.418 217.862 203.123 186.321 185.073 185.495...
                  181.640];
             
%% Theoretical data from the literture
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
ThisWork_HCTH_h2o_h2 = 8.3492.*1e12./Tk.^4 - 91.6939*1e9./Tk.^3 + ...
                  398.4262.*1e6./Tk.^2 - 276.9628.*1e3./Tk + 209.3217;
ThisWork_HCTH_ch4_h2 = 5.0321.*1e12./Tk.^4 - 59.3605*1e9./Tk.^3 + ...
                    270.5051.*1e6./Tk.^2 - 89.2415.*1e3./Tk + -5.1364;
%% Color values
col(1,:)  = [0.655 0.808 0.890]; % Fmd / CO2 / H2O
col(2,:)  = [0.122 0.471 0.706]; % Ftr / CHO-MFR
col(3,:)  = [0.698 0.875 0.541]; % Mch / CHO-H4MPT
col(4,:)  = [0.200 0.627 0.173]; % Mtd / CH-H4MPT
col(5,:)  = [0.984 0.604 0.600]; % Mer / CH2-H4MPT
col(6,:)  = [0.890 0.102 0.110]; % Mtr / CH3-H4MPT
col(7,:)  = [0.992 0.749 0.435]; % Mcr / CH3-S-CoM
col(8,:)  = [1.000 0.498 0.000]; % Mta / CH4
col(9,:)  = [0.792 0.698 0.839]; % Acs / CH3-OH
col(10,:) = [0.416 0.239 0.604]; % Cdh / CH3-COOH
col(11,:) = [0.694 0.349 0.157]; % Frh -- / Hdr -. / CH3-COSCoA
col(12,:) = [0 0 0];             % / F420
col(13,:) = [0.451 0.451 0.451]; % / HS-CoB

%% FIGURE 1 - 2BETA VALUES FOR 1-100 oc
if fig1 == 1
    textLab1 = {'H_2O_{(g)}','CH_4_{(g)}','H_2_{(g)}','H_2O_{(l)}',...
        'CHO-MFR','CHO-H_4MPT','CH-H_4MPT','CH_2-H_4MPT',...
        'CH_3-H_4MPT','CH_3-SCoM','F_{420}H_2','HS-CoB',...
        'CH_3-OH','CH_3-COOH','CH_3-CSCoA'};
    
    x = 1e3./Tk;
    figure(1)
    set(1,'Units','centimeters','Position',[10 5 figWidth figHeight])
%     pos = [1 2 5:8 10:13 15:19];
    pos = [1:8 10:16];
    pos_col = [1,8,2,1,2:7,12:13,9:11];
    for ja = 1:15
        ax1(ja) = plot(x(1:51),betaH(pos(ja),1:51),'Color',col(pos_col(ja),:));
        ax1(ja).LineWidth = figLineWidth;    
        if ja == 3
            ax1(ja).LineStyle = '-.';
        end
        t = text(x(1)+0.01,betaH(pos(ja),1),textLab1{ja});
        t.FontName = 'Helvetica';
        t.FontSize = 12;
        hold on
    end
    box off
%     ax_l = legend(textLab1);
%     ax_l.Location = 'northwest';
    xlabel('1000/T (K)')
    ylabel('^{2}\beta')
    xlim([2.7 3.9])
%     xlim([0 130])
%     ylim([2.5 20])
    grid on
    grid minor
    set(gca,'FontSize',figFontSize)
%     ax_l.FontSize = 11;
    if sav == 1
        print('-f1',[pathstr,'2beta'],'-dpng','-r600');
        print('-f1',[pathstr,'2beta'],'-depsc');
    end
end
%% FIGURE 2A - 2ALPHA EFF VALUES
if fig2 == 1
    textLab2 = {'Ftr';'Mch';'Mtd';'Mer';'Mtr';'Mcr';'Mta';'Acs';'Cdh';...
                'Fmd_p';'Mtd_p';'Mer_p';'Mcr_p'};
    alphaLoc = alphaH(:,1);

    x = betHdat.T_C;
    
    figure(2)
    set(2,'Units','centimeters','Position',[10 5 figWidth figHeight])

    alphalab_leg = {'Fmd','Ftr','Mch','Mtd','Mer',...
        'Mtr','Mcr','Mta','Acs','Cdh','Scheller et al. (2013)'}; 
    pos_s = [2:7,13:15];
    pos_p = [1,9,10,12];
    pos_s_col = 2:10;
    pos_p_col = [1,4,5,7];
    pos_i = [pos_s pos_p];
    pos_col = [pos_s_col pos_p_col];
    
    for i = 1:length(pos_i)
        ax2(i) = plot(x(1:51),1000.*log(alphaH(pos_i(i),1:51)),...
            'Color',col(pos_col(i),:),...
            'LineWidth',figLineWidth);
        hold on
        t = text(103,1000.*log(alphaH(pos_i(i),51)),textLab2{i});
        t.FontName = 'Helvetica';
        t.FontSize = 16;
    end
%     ax2(1) = plot(x(1:51),1000.*log(alphaH(1,1:51)),'Color',col(1,:),...
%         'LineStyle','-.','LineWidth',figLineWidth);
%     hold on
%     for i = 1:length(pos_s)
%         ax2(i+1) = plot(x(1:51),1000.*log(alphaH(pos_s(i),1:51)),'Color',col(pos_s_col(i),:),...
%             'LineWidth',figLineWidth);
%         hold on
%     end
%     for i = 2:length(pos_p)
%         plot(x(1:51),1000.*log(alphaH(pos_p(i),1:51)),'Color',col(pos_p_col(i),:),...
%             'LineStyle','-.','LineWidth',figLineWidth);
%         hold on
%     end
    hold on
    h1 = errorbar(60,16.95,42.67,'o');
    h1.CapSize = 10;
    h1.Color = 'k';
    h1.MarkerFaceColor = col(7,:);
    h1.LineWidth = 1.5;
    h1.MarkerSize = 12;
%     text(63,16.95,'Schller et al. 2013','FontSize',14)
    
%     hl(1).leg = legendflex([ax2 h1], alphalab_leg, 'anchor', [2 2], ...
%         'buffer', [3 -3], ...
%         'ncol', 3, ...
%         'fontsize', 14, ...
%         'xscale', 0.7, ...
%         'box', 'on',...
%         'FontName','Helvetica');
    box off
    xlabel('T (\circC)')
    ylabel(['1000ln^{2}\alpha^{eq} (' char(8240) ')'])
    xlim([0 122])
    ylim([-125 125])
    grid on
    grid minor
    set(gca,'FontSize',figFontSize)
    
    l = legend(h1,'Scheller et al. (2013)');
    l.Location = 'southeast';
    l.Box = 'off';
    if sav == 1
        print('-f2',[pathstr,'2alpha'],'-dpng','-r600');
        print('-f2',[pathstr,'2alpha'],'-depsc');
    end
end

%% FIGURE 2B - 2ALPHA EFF VALUES
if fig9 == 1
    textLab1 = {'CH_4(g)','CHO-MFR','CHO-H_4MPT',...
        'CH-H_4MPT','CH_2-H_4MPT (R)','CH_2-H_4MPT (S)',...
        'CH_3-H_4MPT','CH_3-SCoM','F_{420}H_2','HS-CoB',...
        'CH_3-OH','CH_3-COOH','CH_3-CSCoA'};
    x = betHdat.T_C;

    pos_col = [8,2:5,5,6:7,12:13,9:11];
    pos_i   = [1,3:10,12:14];
    figure(9)
    set(9,'Units','centimeters','Position',[10 5 figWidth figHeight])

    for i = 1:length(pos_i)
        ax2(i) = plot(x(1:51),1000.*log(alphaH_H2O(pos_i(i),1:51)),'Color',col(pos_col(i),:),...
            'LineWidth',figLineWidth);
        hold on
        t = text(103,1000.*log(alphaH_H2O(pos_i(i),51)),textLab1{i});
        t.FontName = 'Helvetica';
        t.FontSize = 13;
    end
    hold on
    
%     hl(1).leg = legendflex([ax2 h1], alphalab_leg, 'anchor', [2 2], ...
%         'buffer', [3 -3], ...
%         'ncol', 3, ...
%         'fontsize', 14, ...
%         'xscale', 0.7, ...
%         'box', 'on',...
%         'FontName','Helvetica');
    box off
    xlabel('T (\circC)')
    ylabel(['1000ln^{2}\alpha^{eq}_{X' char(8211) 'H_2O} (' char(8240) ')'])
    ylim([-230 0])
    xlim([0 150])
    grid on
    grid minor
    set(gca,'FontSize',figFontSize)
    
%     l = legend(h1,'Schller et al. (2013)');
%     l.Location = 'southeast';
%     l.Box = 'off';
    if sav == 1
        print('-f9',[pathstr,'2alpha_h2o'],'-dpng','-r600');
        print('-f9',[pathstr,'2alpha_h2o'],'-depsc');
    end
end

%% FIGURE 3 - 2ALPHA secondary EFF VALUES
if fig3 == 1
    x = 1000./Tk;
    x1a = 1e3./(0+273.15);
    textLab3 = {'Fmd','Mtd','Mer','Mcr','Frh','Mvh/Hdr'};
    col = {[1 1 1].*0;
        [1 1 1].*0;
        [1 1 1].*0;
        [1 1 1].*0;
        [1 1 1].*0;
        [1 1 1].*0};
    sty1 = {'-';
        '-';
        '-';
        '-';
        '-';
        '-'};
    
    figure(3)
    set(3,'Units','centimeters','Position',[10 5 9 9])
    pos = [1 9 10 12 8 11];
    alphaLoc3 = alphaH(:,1);
    alphaLoc3(1) = alphaLoc3(1)*0.95;
    alphaLoc3(5) = alphaLoc3(5)*1.05;
    for ja = 1:6
        plot(x(1:51),1000.*log(alphaH(pos((ja)),1:51)),'Color',col{ja},'LineStyle',sty1{ja})
        %     text(1e3./273.15+0.1,betaLoc(ja,1),textLab1{ja},'FontName','Myriad pro')
        t = text(x1a+0.02,1000.*log(alphaLoc3(pos(ja),1)),textLab3{ja});
        t.FontName = 'Helvetica';
        %     t.Interpreter = 'latex';
        hold on
    end
    box off
    xlabel('1000/T (K)')
    ylabel(['1000ln^{2}\alpha (' char(8240) ')'])
    xlim([2.6 4.05])
    ylim([-850 1000])
    grid on
    grid minor
    if sav == 1
        print('-f3',[pathstr,'2alpha_p'],'-depsc');
    end
end
%% FIGURE 4 - H2O/H2 COMPARISON
if fig4 == 1
    ax1 = figure(4);
    set(4,'Units','centimeters','Position',[10 5 figWidth figHeight])
    % h2o(l)-h2(g)
%     h1 = plot(1000./Tk,Iron19klna_h2o_h2+Horita94_h2o_l_g_klna,'k');
%     h1 = plot(1000./Tk,1000*log(alphaH(15,:)')+Horita94_h2o_l_g_klna,'k');
    h1 = plot(1000./Tk,1000*log(alphaH(19,:)),'k','LineWidth',figLineWidth);
    hold on
    
    h2 = plot(1000./Tk,Suess49klna+Horita94_h2o_l_g_klna,'LineWidth',figLineWidth);
    h2.Color = [227 26 28]./255;%[0.984 0.604 0.600];
    h2.LineStyle = ':';
    hold on
    
    h5 = plot(1000./Tk,Bardo76Calcklna+Horita94_h2o_l_g_klna,'LineWidth',figLineWidth);
    h5.Color = [31 120 180]./255;%[1.000 0.498 0.000];
    h5.LineStyle = '--';
    hold on

    h1a = plot(1000./Tk,ThisWork_HCTH_h2o_h2,':k','LineWidth',figLineWidth);
    hold on
    
    h3 = plot(1000./(Rolston76Tc+273.15),Rolston76klna,'o');
    h3.MarkerSize = marksize;
    h3.MarkerEdgeColor = 'k';
    h3.MarkerFaceColor = 'w';%[0.984 0.604 0.600];
    h3.LineWidth = 0.5;
    hold on
    
    h4 = plot(1000./(Cerrai54Tc(Cerrai54Tc<374)+273.15),Cerrai54klna_l_g(Cerrai54Tc<374),'o');
    h4.MarkerSize = marksize;
    h4.MarkerEdgeColor = 'k';
    h4.MarkerFaceColor = [1 1 1]-0.5;%[1.000 0.498 0.000];
    h4.LineWidth = 0.5;
    hold on
    
    l = legend([h3 h4 h2 h5 h1 h1a],...
        {'Rolston et al. (1976)',...
        'Cerrai et al. (1954)',...
        'Suess (1949)',...
        'Bardo & Wolfsberg (1976)',...
        'This work (M06-L)',...
        'This work (HCTH)'});
    l.Location = 'northwest';
    l.FontSize = 14;
%     l.Box = 'off';
    grid on
%     grid minor
    xlabel('1000/T (K)')
    ylabel(['1000ln^{2}\alpha_{H_2O(l)-H_2(g)} (' char(8240) ')'])
    xlim([1.4 3.8])
    ylim([350 1670])
    set(gca,'FontSize',figFontSize)
    if sav == 1
        %         print('-f4',sprintf('%s2alpha_h2o_h2_comp',pathstr),'-dsvg');
        print('-f4',[pathstr,'2alpha_h2o_h2_comp'],'-dpng','-r600');
        print('-f4',[pathstr,'2alpha_h2o_h2_comp'],'-depsc');
    end
end
%% FIGURE 5 - CH4/H2 COMPARISON
if fig5 == 1
    figure(5)
    set(5,'Units','centimeters','Position',[10 5 figWidth figHeight])
    
    h6 = plot(1e3./Tk,Bottinga69aklna,'LineWidth',figLineWidth);
    h6.Color = [227 26 28]./255;%[0.694 0.349 0.157];
    h6.LineStyle = '-';
    hold on
    
    h7 = plot(1e3./Tk,Richet77aklna,'LineWidth',figLineWidth);
    h7.Color = [0.122 0.471 0.706];
    h7.LineStyle = '-';
    hold on
%     
%     h8 = plot(1e3./Tk,Richet77bklna);
%     h8.Color = [0.122 0.471 0.706];
%     h8.LineStyle = '-';
%     hold on
    
    h9 = plot(1e3./Tk,Horibe95Calcklna,'LineWidth',figLineWidth);
    h9.Color = [51 160 44]./255;%[0.416 0.239 0.604];
    h9.LineStyle = '-';
    hold on

    h5 = plot(1e3./Tk,1000*log(alphaH(17,:)),'LineWidth',figLineWidth);
    h5.Color = 'k';
    hold on
    
    h5a = plot(1e3./Tk,ThisWork_HCTH_ch4_h2,'LineWidth',figLineWidth);
    h5a.Color = 'k';
    h5a.LineStyle = ':';
    hold on
    
    h10 = plot(1e3./(Horibe1995Tc+273.15),Horibe1995klna,'o');
    h10.MarkerEdgeColor = 'k';
    h10.MarkerFaceColor = 'w';%[51 160 44]./255;%[0.416 0.239 0.604];
    h10.MarkerSize = marksize;
    h10.LineWidth = 0.5;
    hold on

%     turner_xdat = [2.115 2.203 2.33 2.496 2.683 2.866 3.085 3.243 3.356 3.476 3.533 3.617];
%     turner_ydat = [570.874 602.589 647.896 724.919 826.861 906.149 1028.479 1098.706 1153.074 1209.709 1255.016 1304.854];
%     h11 = plot(turner_xdat,turner_ydat,'o');
%     h11.MarkerEdgeColor = 'k';
%     h11.MarkerFaceColor = [0.694 0.349 0.157];
%     h11.MarkerSize = marksize;
%     h11.LineWidth = 0.5;
%     hold on
    
    l = legend([h10 h9 h6 h7 h5 h5a],...
        {'Horibe & Craig (1995)',...
        'Horibe & Craig (1995)',...
        'Bottinga (1969)',...
        'Richet et al. (1977)',...
        'This work (M06-L)',...
        'This work (HCTH)'});
    l.Location = 'northwest';
%     l.Box = 'off';
    l.FontSize = 14;
    % l.EdgeColor = 'w';
    grid on
%     grid minor
    xlabel('1000/T (K)')
    ylabel(['1000ln^{2}\alpha_{CH_4(g)-H_2(g)} (' char(8240) ')'])
    xlim([1 3.8])
    
%     xti = [0.5 2.5 5 7.5 10 12.5];
%     set(gca,'XTick',xti)
%     set(gca,'XTickLabel',round(sqrt(1e6./xti)-273.15))
    set(gca,'FontSize',figFontSize)
    if sav == 1
%         print('-f5',sprintf('%s2alpha_ch4_h2_comp',pathstr),'-dpng','-r600');
%         print('-f5',sprintf('%s2alpha_ch4_h2_comp',pathstr),'-dsvg');
        print('-f5',[pathstr,'2alpha_ch4_h2_comp'],'-dpng','-r600');
        print('-f5',[pathstr,'2alpha_ch4_h2_comp'],'-depsc');
    end
end
%% FIGURE 6 - CH4/H2O COMPARISON
if fig6 == 1
    figure(6)
    set(6,'Units','centimeters','Position',[10 5 figWidth figHeight])
    % h2o(l)-h2(g)
%     h1 = plot(Tc(Tc<=374),-Horita94_h2o_l_g_klna(Tc<=374)'+1000.*log(alphaH(13,(Tc<=374))),'k');
    h1 = plot(Tc,1000.*log(alphaH(16,:)),'k','LineWidth',figLineWidth);
    hold on
    h2 = plot(Tc(Tc<=374),-Horita94_h2o_l_g_klna(Tc<=374)-Bottinga69bklna(Tc<=374),...
        'Color',[0.694 0.349 0.157],'LineWidth',figLineWidth);
    h3 = plot(Tc(Tc<=374),-Horita94_h2o_l_g_klna(Tc<=374)-Suess49klna(Tc<=374)+...
        Horibe95Calcklna(Tc<=374),'LineWidth',figLineWidth);
    h3.Color = [0.890 0.102 0.110];
    h3.LineStyle = '-';
    hold on
    h4 = plot(Tc(Tc<=374),-Horita94_h2o_l_g_klna(Tc<=374)-Cerrai54Calcklna(Tc<=374)+...
        Horibe95Calcklna(Tc<=374),'LineWidth',figLineWidth);
    h4.Color = [1.000 0.498 0.000];
    h4.LineStyle = '-';
    h5 = plot(Tc(Tc<=200),-Rolston76Calcklna(Tc<=200)+Horibe95Calcklna(Tc<=200),...
        'LineWidth',figLineWidth);
    h5.Color = [0.984 0.604 0.600];
    h5.LineStyle = '-';
    h6 = plot(Tc(Tc<=374),-Horita94_h2o_l_g_klna(Tc<=374)-Richet77cklna(Tc<=374),...
        'Color',[0.122 0.471 0.706],'LineWidth',figLineWidth);
    
    l = legend([h3 h4 h5 h2 h6 h1],...
        {'HW94+S49+HC95',...
        'HW94+C54+HC95',...
        'R76+HC95',...
        'HW94+B69',...
        'HW94+R77',...
        'HW94+G20'});
    l.Location = 'southeast';
    l.Box = 'on';
    l.FontSize = 14;
    grid on
%     grid minor
    xlabel('T (\circC)')
    ylabel(['1000ln^2\alpha_{CH_4-H_2O(l)} (' char(8240) ')'])
    xlim([-10 384])
    ylim([-310 -60])
    set(gca,'FontSize',figFontSize)
    if sav == 1
%         print('-f6',sprintf('%s2alpha_h2o_ch4_comp',pathstr),'-dpng','-r600');
%         print('-f6',sprintf('%s2alpha_h2o_ch4_comp',pathstr),'-dsvg');
        print('-f6',[pathstr,'2alpha_h2o_ch4_comp'],'-dpng','-r600');
        print('-f6',[pathstr,'2alpha_h2o_ch4_comp'],'-depsc');
    end
end

%% FIGURE 7 - H2O/H2 COMPARISON - fully substituted
if fig7 == 1
    figure(7)
    set(7,'Units','centimeters','Position',[10 5 figWidth figHeight])
    
    local_path = '/Users/jonathag/Dropbox (Weizmann Institute)/Jonathan_with_Itay/Paper 2 EFFs/';
    file_name  = 'alpha-beta_fully_substituted.xlsx';
    betHdat_full_subs    = readtable([local_path,file_name],'Sheet',2);
    
    h3 = plot(1000./(Rolston76Tc+273.15),Rolston76klna,'o');
    h3.MarkerSize = marksize;
    h3.MarkerEdgeColor = 'k';
    h3.MarkerFaceColor = [0.984 0.604 0.600];
    h3.LineWidth = 0.5;
    hold on
    
    h4 = plot(1000./(Cerrai54Tc(Cerrai54Tc<374)+273.15),Cerrai54klna_l_g(Cerrai54Tc<374),'o');
    h4.MarkerSize = marksize;
    h4.MarkerEdgeColor = 'k';
    h4.MarkerFaceColor = [1.000 0.498 0.000];
    h4.LineWidth = 0.5;
    hold on

    h1 = plot(1000./Tk,1000*log(betHdat_full_subs.HDO_l./...
        betHdat_full_subs.HD_g),':k','LineWidth',figLineWidth);
    hold on
    
    h2 = plot(1000./Tk,1000*log(alphaH(19,:)),'k','LineWidth',figLineWidth);
    hold on
        
    l = legend([h3 h4 h2 h1],...
        {'Rolston et al., 1976',...
        'Cerrai et al., 1954',...
        'Singly substituted',...
        'Fully substituted'});
    l.Location = 'northwest';
    l.FontSize = 14;
%     l.Box = 'off';
    grid on
%     grid minor
    xlabel('1000/T (K)')
    ylabel(['1000ln^{2}\alpha_{H_2O(l)-H_2(g)} (' char(8240) ')'])
    xlim([1.4 3.8])
    set(gca,'FontSize',figFontSize)
    if sav == 1
        %         print('-f4',sprintf('%s2alpha_h2o_h2_comp',pathstr),'-dsvg');
        print('-f7',[pathstr,'sup_2alpha_h2o_h2_comp'],'-dpng','-r600');
        print('-f7',[pathstr,'sup_2alpha_h2o_h2_comp'],'-depsc');
    end
end

%% FIGURE 8 - CH4/H2 COMPARISON fully substituted
if fig8 == 1
    figure(8)
    set(8,'Units','centimeters','Position',[10 5 figWidth figHeight])
    
    local_path = '/Users/jonathag/Dropbox (Weizmann Institute)/Jonathan_(with_Itay)/Paper 2 EFFs/';
    file_name  = 'alpha-beta_fully_substituted.xlsx';
    betHdat_full_subs    = readtable([local_path,file_name],'Sheet',2);

    h4 = plot(1e3./(Horibe1995Tc+273.15),Horibe1995klna,'o');
    h4.MarkerEdgeColor = 'k';
    h4.MarkerFaceColor = [0.416 0.239 0.604];
    h4.MarkerSize = marksize;
    h4.LineWidth = 0.5;
    hold on

    turner_xdat = [2.115 2.203 2.33 2.496 2.683 2.866 3.085 3.243 3.356 3.476 3.533 3.617];
    turner_ydat = [570.874 602.589 647.896 724.919 826.861 906.149 1028.479 1098.706 1153.074 1209.709 1255.016 1304.854];
    h5 = plot(turner_xdat,turner_ydat,'o');
    h5.MarkerEdgeColor = 'k';
    h5.MarkerFaceColor = [0.694 0.349 0.157];
    h5.MarkerSize = marksize;
    h5.LineWidth = 0.5;
    hold on
    
    h2 = plot(1e3./Tk,1000*log(alphaH(17,:)),'LineWidth',figLineWidth);
    h2.Color = 'k';
    hold on
    
    h3 = plot(1000./Tk,1000*log(betHdat_full_subs.C12H3D_g./...
        betHdat_full_subs.HD_g),':k','LineWidth',figLineWidth);
    
    l = legend([h5 h4 h2 h3],...
        {'Turner et al., 2020',...
        'Horibe & Craig, 1995',...
        'Singly substituted',...
        'Fully substituted'});
    l.Location = 'northwest';
    l.FontSize = 14;
    grid on
    xlabel('1000/T (K)')
    ylabel(['1000ln^{2}\alpha_{CH_4(g)-H_2(g)} (' char(8240) ')'])
    xlim([1 3.8])
    
    set(gca,'FontSize',figFontSize)
    if sav == 1
        print('-f8',[pathstr,'sup_2alpha_ch4_h2_comp'],'-dpng','-r600');
        print('-f8',[pathstr,'sup_2alpha_ch4_h2_comp'],'-depsc');
    end
end
end
