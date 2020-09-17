% 19.11.18
% EFF paper
% Plot beta values for hydrogen isotopes

clear
% [betH,betHlab] = xlsread('C:\Jonathan\Dropbox (Weizmann Institute)\Paper 2 Methanogenesis EFFs\alpha-beta.xlsx',4);
% betH(1,:) = [];

local_path = '/Users/jonathag/Dropbox (Weizmann Institute)/Paper 2 Methanogenesis EFFs/';
pathstr    = [local_path,'Figures/']; 
file_name  = 'alpha-beta.xlsx';
betHdat    = readtable([local_path,file_name],'Sheet',2);

fig1 = 0; % 2beta
fig2 = 0; % 1000ln2a secondary
fig3 = 0; % 1000ln2a primary
fig4 = 0; % H2O(l)-H2 comparison
fig5 = 1; % CH4-H2 comparison
fig6 = 0; % CH4-H2O comparison
sav  = 0;
plot_EFF_Hfig(betHdat,fig1,fig2,fig3,fig4,fig5,fig6,sav,pathstr)
% close all

function [] = plot_EFF_Hfig(betHdat,fig1,fig2,fig3,fig4,fig5,fig6,sav,pathstr)
marksize = 5;
addpath(genpath('/Users/jonathag/Dropbox (Weizmann Institute)/Jonathan_(with_Itay)/Code/legendflex-pkg-master'))

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
% Horita94_h2o_l_g_klna(Tc>374) = nan;
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
col(11,:) = [0.694 0.349 0.157]; % Frh -- / Hdr -. / CH3-COSCoA

%% FIGURE 1 - 2BETA VALUES FOR 1-100 oc
if fig1 == 1
    textLab1 = {'H_2O(g)','CH_4(g)','H_2(g)','CHO-MFR','CHO-H_4MPT',...
        'CH-H_4MPT(R)','CH-H_4MPT(S)','CH_2-H_4MPT',...
        'CH_3-H_4MPT','CH_3-S-CoM','F_{420}H_2','HS-CoB',...
        'Methanol','Acetate','Acetyl-CoA'};
    textLab1 = {'H2O(g)','CH4(g)','H2(g)','CHO-MFR','CHO-H4MPT',...
        'CH-H4MPT(R)','CH-H4MPT(S)','CH2-H4MPT',...
        'CH3-H4MPT','CH3-S-CoM','F420H2','HS-CoB',...
        'CH3-OH','CH3-COOH','CH3-CSCoA'};
    betaLoc = betaH(:,1);
    
    x1 = 1e3./Tk;
    x1a = 1e3./(0+273.15);
    x2a = 1e3./(100+273.15);
    
    figure(1)
    set(1,'Units','centimeters','Position',[10 5 9 9])
%     pos = [1 2 5:8 10:13 15:19];
%     pos = [1:3 5:16];
    pos_col = [1,8,2,2:5,5,6:7,3:4,9:11];
    for ja = 1:15
        ax1(ja) = plot(x1(1:51),betaH((ja),1:51),'Color',col(pos_col(ja),:));
%         t = text(x1a+0.02,betaLoc(pos(ja),1),textLab1{ja});
%         t.FontName = 'Helvetica';
%         t.FontSize = 8;
        if ja == 1 || ja == 2 || ja == 3 || ja == 11 || ja == 12
            ax1(ja).LineStyle = '-.';
        end
        hold on
    end
    box off
%     legendflex(ax1, textLab1, 'anchor', [3 3], ...
%         'buffer', [-3 -3], ...
%         'ncol', 3, ...
%         'fontsize', 8, ...
%         'xscale', 0.7, ...
%         'box', 'on',...
%         'FontName','Helvetica')
    ax_l = legend(textLab1);
    ax_l.Location = 'eastoutside';
    xlabel('1000/T (K)')
    ylabel('^{2}\beta')
%     xlim([2.6 4.2])
%     xlim([0 100])
%     ylim([14 19])
    grid on
    grid minor
    if sav == 1
        print('-f1',[pathstr,'2beta'],'-dpng','-r600');
        print('-f1',[pathstr,'2beta'],'-depsc');
    end
end
%% FIGURE 2 - 2ALPHA secondary EFF VALUES
if fig2 == 1
    textLab2 = {'Fmd';'Ftr';'Mch';'Mtd';'Mer';'Mtr';'Mcr';'Mta';'Acs';'Cdh'};
    alphaLoc = alphaH(:,1);

%     x1 = 1e3./Tk;
%     x1a = 1e3./(0+273.15);
    x1 = betHdat.T_C;
    
    figure(2)
    set(2,'Units','centimeters','Position',[10 5 9 9])

    alphalab_leg = {'Fmd','Ftr','Mch','Mtd','Mer','Mtr','Mcr','Mta','Acs','Cdh'}; 
    pos_s = [2:7,13:15];
    pos_p = [1,9,10,12];
    pos_s_col = 2:10;
    pos_p_col = [1,4,5,7];
    
    ax2(1) = plot(x1(1:51),1000.*log(alphaH(1,1:51)),'Color',col(1,:),'LineStyle','-.');
    hold on
    for i = 1:length(pos_s)
        ax2(i+1) = plot(x1(1:51),1000.*log(alphaH(pos_s(i),1:51)),'Color',col(pos_s_col(i),:));
        hold on
    end
    for i = 2:length(pos_p)
        plot(x1(1:51),1000.*log(alphaH(pos_p(i),1:51)),...
            'Color',col(pos_p_col(i),:),'LineStyle','-.');
        hold on
    end
    hold on
    h1 = errorbar(60,16.95,42.67,'o');
    h1.CapSize = 1;
    h1.Color = 'k';
    h1.MarkerFaceColor = [1 1 1].*0.9;
    
    hl(1).leg = legendflex(ax2, alphalab_leg, 'anchor', [2 2], ...
        'buffer', [3 -3], ...
        'ncol', 3, ...
        'fontsize', 8, ...
        'xscale', 0.7, ...
        'box', 'on',...
        'FontName','Helvetica');
    box off
    xlabel('T (\circC)')
    ylabel(['1000ln^{2}\alpha (' char(8240) ')'])
    ylim([-140 210])
    grid on
    grid minor
        
    l = legend(h1,'Schller et al. (2013)');
    l.Location = 'southeast';
    l.Box = 'off';
    if sav == 1
        print('-f2',[pathstr,'2alpha'],'-dpng','-r600');
        print('-f2',[pathstr,'2alpha'],'-depsc');
    end
end
%% FIGURE 3 - 2ALPHA secondary EFF VALUES
if fig3 == 1
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
        plot(x1(1:51),1000.*log(alphaH(pos((ja)),1:51)),'Color',col{ja},'LineStyle',sty1{ja})
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
    figure(4)
    set(4,'Units','centimeters','Position',[10 5 9 9])
    % h2o(l)-h2(g)
%     h1 = plot(1000./Tk,Iron19klna_h2o_h2+Horita94_h2o_l_g_klna,'k');
%     h1 = plot(1000./Tk,1000*log(alphaH(15,:)')+Horita94_h2o_l_g_klna,'k');
    h1 = plot(1000./Tk,1000*log(alphaH(19,:)),'k');
    hold on
    h2 = plot(1000./Tk,Suess49klna+Horita94_h2o_l_g_klna);
    h2.Color = [0.890 0.102 0.110];
    h2.LineStyle = '-';
    hold on
    h3 = plot(1000./(Rolston76Tc+273.15),Rolston76klna,'o');
    h3.MarkerSize = marksize;
    h3.MarkerEdgeColor = 'k';
    h3.MarkerFaceColor = [0.984 0.604 0.600];
    h3.LineWidth = 0.5;
    hold on
    h4 = plot(1000./(Cerrai54Tc(Cerrai54Tc<700)+273.15),Cerrai54klna_l_g(Cerrai54Tc<700),'o');
    h4.MarkerSize = marksize;
    h4.MarkerEdgeColor = 'k';
    h4.MarkerFaceColor = [1.000 0.498 0.000];
    h4.LineWidth = 0.5;
    hold on
    
    l = legend([h4 h3 h2 h1],...
        {'Cerrai et al., 1954',...
        'Rolston et al., 1976',...
        'Suess, 1949',...
        'This work'});
    l.Location = 'northwest';
    l.Box = 'off';
    grid on
    grid minor
    xlabel('1000/T (K)')
    ylabel(['1000ln^{2}\alpha_{H_2O(l)-H_2(g)} (' char(8240) ')'])
    if sav == 1
        %         print('-f4',sprintf('%s2alpha_h2o_h2_comp',pathstr),'-dsvg');
        print('-f4',[pathstr,'2alpha_h2o_h2_comp'],'-dpng','-r600');
        print('-f4',[pathstr,'2alpha_h2o_h2_comp'],'-depsc');
    end
end
%% FIGURE 5 - CH4/H2 COMPARISON
if fig5 == 1
    figure(5)
    set(5,'Units','centimeters','Position',[10 5 9 9])
    
    h6 = plot(1e3./Tk,Bottinga69aklna);
    h6.Color = [0.694 0.349 0.157];
    h6.LineStyle = '-';
    hold on
    
    h7 = plot(1e3./Tk,Richet77aklna);
    h7.Color = [0.122 0.471 0.706];
    h7.LineStyle = '-';
    hold on
%     
%     h8 = plot(1e3./Tk,Richet77bklna);
%     h8.Color = [0.122 0.471 0.706];
%     h8.LineStyle = '-';
%     hold on
    
    h9 = plot(1e3./Tk,Horibe95Calcklna);
    h9.Color = [0.416 0.239 0.604];
    h9.LineStyle = '-';
    hold on

    h5 = plot(1e3./Tk,1000*log(alphaH(17,:)));
    h5.Color = 'k';
    hold on
    
    h10 = plot(1e3./(Horibe1995Tc+273.15),Horibe1995klna,'o');
    h10.MarkerEdgeColor = 'k';
    h10.MarkerFaceColor = [0.416 0.239 0.604];
    h10.MarkerSize = marksize;
    h10.LineWidth = 0.5;
    hold on

    turner_xdat = [2.115 2.203 2.33 2.496 2.683 2.866 3.085 3.243 3.356 3.476 3.533 3.617];
    turner_ydat = [570.874 602.589 647.896 724.919 826.861 906.149 1028.479 1098.706 1153.074 1209.709 1255.016 1304.854];
    h11 = plot(turner_xdat,turner_ydat,'o');
    h11.MarkerEdgeColor = 'k';
    h11.MarkerFaceColor = [0.694 0.349 0.157];
    h11.MarkerSize = marksize;
    h11.LineWidth = 0.5;
    hold on
    
    l = legend([h10 h11 h9 h6 h7 h5],...
        {'Horibe & Craig, 1995',...
        'Turner et al., 2020',...
        'Horibe & Craig, 1995',...
        'Bottinga, 1969',...
        'Richet et al., 1977',...
        'This work'});
    l.Location = 'northwest';
    l.Box = 'off';
    l.FontSize = 8;
    % l.EdgeColor = 'w';
    grid on
    grid minor
    xlabel('1000/T (K)')
    ylabel(['1000ln^{2}\alpha_{CH_4(g)-H_2(g)} (' char(8240) ')'])
    xlim([0.8 3.8])
    
%     xti = [0.5 2.5 5 7.5 10 12.5];
%     set(gca,'XTick',xti)
    % set(gca,'XTickLabel',round(sqrt(1e6./xti)-273.15))
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
    set(6,'Units','centimeters','Position',[10 5 9 9])
    % h2o(l)-h2(g)
%     h1 = plot(Tc(Tc<=374),-Horita94_h2o_l_g_klna(Tc<=374)'+1000.*log(alphaH(13,(Tc<=374))),'k');
    h1 = plot(Tc,1000.*log(alphaH(16,:)),'k');
    hold on
    h2 = plot(Tc(Tc<=374),-Horita94_h2o_l_g_klna(Tc<=374)-Bottinga69bklna(Tc<=374),'Color',[0.694 0.349 0.157]);
    h3 = plot(Tc(Tc<=374),-Horita94_h2o_l_g_klna(Tc<=374)-Suess49klna(Tc<=374)+Horibe95Calcklna(Tc<=374));
    h3.Color = [0.890 0.102 0.110];
    h3.LineStyle = '-';
    hold on
    h4 = plot(Tc(Tc<=374),-Horita94_h2o_l_g_klna(Tc<=374)-Cerrai54Calcklna(Tc<=374)+Horibe95Calcklna(Tc<=374));
    h4.Color = [1.000 0.498 0.000];
    h4.LineStyle = '-';
    h5 = plot(Tc(Tc<=200),-Rolston76Calcklna(Tc<=200)+Horibe95Calcklna(Tc<=200));
    h5.Color = [0.984 0.604 0.600];
    h5.LineStyle = '-';
    h6 = plot(Tc(Tc<=374),-Horita94_h2o_l_g_klna(Tc<=374)-Richet77cklna(Tc<=374),'Color',[0.122 0.471 0.706]);
    
    l = legend([h3 h4 h5 h2 h6 h1],...
        {'H94+S49+H95',...
        'H94+C54+H95',...
        'R76+H95',...
        'H94+B69',...
        'H94+R77',...
        'H94+G20'});
    l.Location = 'southeast';
    l.Box = 'on';
    grid on
    grid minor
    xlabel('T (\circC)')
    ylabel(['1000ln^2\alpha_{CH_4-H_2O(l)} (' char(8240) ')'])
    xlim([-10 384])
    ylim([-310 -60])
    if sav == 1
%         print('-f6',sprintf('%s2alpha_h2o_ch4_comp',pathstr),'-dpng','-r600');
%         print('-f6',sprintf('%s2alpha_h2o_ch4_comp',pathstr),'-dsvg');
        print('-f6',[pathstr,'2alpha_h2o_ch4_comp'],'-dpng','-r600');
        print('-f6',[pathstr,'2alpha_h2o_ch4_comp'],'-depsc');
    end
end

end
