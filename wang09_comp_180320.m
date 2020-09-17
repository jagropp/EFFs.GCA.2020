% 20.3.18
% Compilation of results from the Wang papers

% 1. 2,4-dimethyl-3-petanone 3o Ha
% 2. 3,5-dimethyl-4-heptanone 3o Ha
% 3. 2-methyl-3-hexanone 3o Ha
% 4. 2-methyl-3-hexanone 2o Ha
% 5. 4-heptanone 2o Ha
% 6. 5-nonanone 2o Ha
% 7. 2-heptanone 1o Ha
% 8. 2-heptanone 2o Ha
% 9. cyclohexanone 2o Heql
% 10. cyclohexanone 2o Haxi
% 11. H2O(l)

% 1. 2-methyl-cyclohexanone 2o Ha eql
% 2. 2-methyl-cyclohexanone 2o Ha axl
% 3. 2-methyl-cyclohexanone 3o Ha axl
% 4. 3-methyl-cyclohexanone 2o Ha eql
% 5. 3-methyl-cyclohexanone 2o Ha axl
% 6. 4-methyl-cyclohexanone 2o Ha eql
% 7. 4-methyl-cyclohexanone 2o Ha axl
% 8. Cis-2,6-dimethyl-cyclohexanone 3o Ha axl
% 9. 2,2,6-trimethyl-cyclohexanone 3o Ha axl
% 10. H2O(l)

Tc09 = [0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100];
bet09 = [20.215 18.851 17.625 16.518 15.517 14.608 13.782 13.028 12.338 11.706 11.126 10.592 10.099 9.643 9.221 8.83 8.467 8.128 7.813 7.518 7.243;
    19.655 18.339 17.154 16.086 15.118 14.239 13.44 12.71 12.043 11.431 10.868 10.35 9.872 9.431 9.021 8.642 8.289 7.96 7.654 7.367 7.1;
    20.292 18.921 17.689 16.577 15.571 14.658 13.828 13.07 12.378 11.743 11.16 10.623 10.128 9.671 9.247 8.854 8.489 8.15 7.833 7.537 7.261;
    16.812 15.754 14.797 13.931 13.143 12.425 11.769 11.168 10.617 10.109 9.642 9.21 8.81 8.439 8.094 7.774 7.475 7.196 6.935 6.691 6.462;
    16.817 15.759 14.803 13.936 13.149 12.431 11.775 11.174 10.623 10.115 9.647 9.215 8.815 8.444 8.1 7.779 7.48 7.201 6.94 6.696 6.467;
    16.671 15.625 14.68 13.823 13.043 12.333 11.684 11.09 10.544 10.042 9.579 9.151 8.755 8.387 8.046 7.728 7.432 7.156 6.897 6.655 6.428;
    15.753 14.799 13.934 13.148 12.432 11.777 11.178 10.628 10.122 9.656 9.225 8.826 8.456 8.112 7.792 7.494 7.216 6.955 6.712 6.483 6.268;
    16.819 15.761 14.805 13.938 13.15 12.431 11.775 11.174 10.623 10.115 9.647 9.215 8.815 8.444 8.1 7.779 7.48 7.201 6.94 6.696 6.467;
    18.215 17.038 15.976 15.016 14.144 13.351 12.627 11.966 11.359 10.802 10.289 9.815 9.378 8.973 8.597 8.247 7.922 7.618 7.335 7.07 6.821;
    16.881 15.816 14.853 13.98 13.188 12.465 11.805 11.201 10.647 10.137 9.667 9.233 8.831 8.459 8.113 7.791 7.491 7.211 6.949 6.704 6.474;
    19.291 18.053 16.939 15.932 15.02 14.192 13.436 12.747 12.114 11.533 10.998 10.504 10.048 9.625 9.232 8.867 8.527 8.209 7.912 7.634 7.374];
Tc13 = [0 10 20 25 30 40 50 60 70 80 90 100];
bet13 = [18.303 16.05 14.206 13.408 12.68 11.404 10.328 9.412 8.626 7.948 7.358 6.842;
    16.97 14.928 13.252 12.525 11.861 10.695 9.709 8.868 8.145 7.52 6.975 6.497;
    18.692 16.349 14.436 13.61 12.858 11.541 10.433 9.491 8.686 7.991 7.389 6.862;
    18.159 15.928 14.101 13.311 12.59 11.326 10.259 9.351 8.573 7.9 7.315 6.803;
    16.818 14.798 13.14 12.421 11.764 10.611 9.635 8.802 8.087 7.467 6.927 6.454;
    18.159 15.928 14.103 13.312 12.591 11.327 10.261 9.353 8.574 7.902 7.317 6.805;
    16.804 14.787 13.132 12.414 11.758 10.605 9.631 8.799 8.084 7.465 6.926 6.454;
    18.875 16.502 14.566 13.73 12.969 11.637 10.516 9.564 8.75 8.049 7.44 6.908;
    19.082 16.676 14.713 13.866 13.095 11.746 10.611 9.648 8.824 8.114 7.498 6.961;
    19.291 16.939 15.02 14.192 13.436 12.114 10.998 10.048 9.232 8.527 7.912 7.374];

alpha09([1 2 4 5],:) = bet09([1 2 5 6],:)./bet09(end,:);
alpha09(3,:) = (bet09(3,:)+2.*bet09(4,:))./3./bet09(end,:);
alpha09(6,:) = (3.*bet09(6,:)+2.*bet09(7,:))./5./bet09(end,:);
% alpha09(7,:) = ([bet09(8,:);bet09(9,:)])./bet09(end,:);
% alpha09 = bet09(1:end-1,:)./bet09(end,:);
alpha13 = bet13(1:end-1,:)./bet13(end,:);

% Wang 2009
% 1. 2,4-Dimethyl-3-pentanone
% 2. 3,5-Dimethyl-4-heptanone
% 3. 2-Methyl-3-hexanone
% 4. 4-Heptanone
% 5. 5-Nonanone
% 6. 2-Heptanone
% 7. Cyclohexanone
% Wang 2013
% 8. Cyclohexanone
% 9. 2-Methylcyclohexanone
% 10. 3-Methylcyclohexanone
% 11. 4-Methylcyclohexanone
% 12. 2,6-Dimethylcyclohexanone
% 13. 2,2,6-Trimethylcyclohexanone

TcExp = [25 50 70];
epsExp = [27 4 -25; 8 -18 -32; -81 -65 -90; -130 -129 -112; -122 -112 -128;
    -152 -145 -140; -161 -155 -138; -148 -158 -152; -102 -105 -103;
    -148 -154 -139; -133 -142 -137; -8 -30 -42; 7 -31 -34];

errExp = [11 12 20; 9 20 17; 6 8 17; 22 17 21; 12 18 29; 30 25 19; 11 10 10;
          11 7 6; 7 11 6; 5 13 9; 10 10 10; 12 10 18; 26 19 19];

col = [166,206,227;
     31,120,180;
    178,223,138;
     51,160,44;
    251,154,153;
    227, 26, 28;
    253,191,111;
    255,127,  0;
    202,178,214;
    106, 61,154]./255;

% pos = [1 2 3 5 6 7 9];
alkName = {'2,4-Dimethyl-3-pentanone','3,5-Dimethyl-4-heptanone',...
    '2-Methyl-3-hexanone','4-Heptanone','5-Nonanone','2-Heptanone'};

h = figure;
set(h,'units','centimeters','position',[10 10 10 8])
% subplot(1,2,1)    
h1 = plot(Tc09,1000.*(alpha09-1),'--');
set(h1,{'Color'},num2cell(col(1:6,:),2));
hold on
for i = 1:6
    errorbar(TcExp,epsExp(i,:),errExp(i,:),'Color',col(i,:))
    hold on
    scatter(TcExp,epsExp(i,:),'filled','MarkerFaceColor',col(i,:),'MarkerEdgeColor','k')
    text(73,epsExp(i,3),alkName{i},'FontName','Myriad pro')
    hold on
end
xlabel('T [$^{\circ}$C]','interpreter','latex')
ylabel('1000$\times$ln$\alpha$(alkane-H$_2$O$_{(l)}$)','interpreter','latex')
% subplot(1,2,2)
% plot(Tc13,1000.*(alpha13-1))
% hold on
% for i = 8:13
%     scatter(TcExp,epsExp(i,:),'filled')
%     hold on
% end
% xlabel('T [$^{\circ}$C]','interpreter','latex')
% ylabel('1000$\times$ln$\alpha$(alkane-H$_2$O$_{(l)}$)','interpreter','latex')

% M06-L gas phase 25�C
alphaDFT(:,1) = [1.0185 0.9925 0.9361 0.8936 0.8962 0.8629];
% M06-L SMD water 25�C
alphaDFT(:,2) = [1.0329 0.9989 0.9390 0.9082 NaN    0.8645];
% SOGGA11 gas phase
alphaDFT(:,3) = [0.9854 1.0051 0.9008 0.8274 0.8851 0.8055];

hold on
for i = 1:6
    scatter(25,1000.*(alphaDFT(i,3)-1),'filled','d','MarkerFaceColor',col(i,:),'MarkerEdgeColor','k');
    hold on
end
% scatter(repmat(25,length(1000.*(alphaDFT-1)),1),1000.*(alphaDFT(:,2)-1),'filled');    











