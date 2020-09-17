% 12.3.18
% Data from Mark compared to Wang 2009, based on:
% 180221 Equilibrium Isotopic Fractionation - Liu 2.

% Molecules (each line in aeq):
% 2,2,6-trimethylcyclohexanone
% 2,4-dimethyl3pentanone
% 2-heptanone
% 2-methyl3hexanone
% 2-methylcyclohexanone
% 3,5-dimethyl4heptanone
% 3-methylcyclohexanone
% 4-heptanone
% 4-methylcyclohexanone
% 9-nonanone (5-nonanone)
% cyclohexanone

% Methods (each column in aeq)
% 1. Wang (exp)
% 2. Wang (calc) 
% 3. M06-L (gp) 
% 4. M06-L (SMD) 
% 5. SOGGA11 (gp) 
% 6. OLYP (gp) 
% 7. M11-L (gp)
% 8. Experimental uncertainties

aeq = ...
[1.007 1.054 0.974 0.996 0.918 0.968 0.837 0.026;
 1.027 1.082 1.019 1.044 0.985 1.027 0.879 0.011;
 0.848 0.915 0.863 0.869 0.806 0.856 0.747 0.030;
 0.919 1.001 0.936 0.946 0.901 0.935 0.806 0.006;
 0.898 1.002 0.945 0.949 0.895 0.935 0.814 0.007;
 1.008 1.082 0.992 1.009 1.005 1.033 0.857 0.016;
 0.852 0.978 0.929 0.929 0.889 0.919 0.798 0.005;
 0.870 0.945 0.894 0.915 0.827 0.867 0.776 0.022;
 0.867 0.978 0.931 0.929 0.886 0.919 0.803 0.010;
 0.878 0.937 0.896 NaN   0.885 0.907 0.775 0.012;
 0.852 0.981 0.932 0.930 0.884 0.919 0.802 0.011];

aeq(:,2) = [0.977;
1.003;
0.848;
0.928;
0.929;
1.003;
0.907;
0.876;
0.906;
0.869;
0.910]; % a manual correction I added, I think there's an error in Mark's 
        % calculations for the beta of water he uses.

h1 = figure;
set(h1,'units','centimeters','position',[20 10 10 8])
for i = 2:7
    scatter(aeq(:,1),aeq(:,i),'filled')
    hold on
end
hold on
% line([min(min(aeq))-0.05 max(max(aeq))+0.05],[min(min(aeq))-0.05 max(max(aeq))+0.05])
% xlim([min(min(aeq(:,1)))-0.05 max(max(aeq(:,1)))+0.05])
% ylim([min(min(aeq))-0.05 max(max(aeq))+0.05])
xlim([0.7 1.1])
ylim([0.7 1.1])
line([0.7 1.1],[0.7 1.1])

h2 = figure;
set(h2,'units','centimeters','position',[20 10 10 8])
x = repmat(linspace(1,11,11),7,1);
for i = 1:7
    scatter(x(i,:),aeq(:,i),'filled')
%     plot(x(i,:),aeq(:,i))
    hold on
end
errorbar(linspace(1,11,11),aeq(:,1),aeq(:,8),'LineStyle','none')
xlim([0 12])
ylabel('1000$\times$ln$\alpha$(Org-H$_2$O$_{(\mathrm{liq})}$)','interpreter','latex')
legend('Wang 09 Exp','Wang 09 B3LYP','M06-L (gp)','M06-L (SMD)','SOGGA11 (gp)','OLYP (gp)','M11-L (gp)','Location','best')
legend('boxoff')




