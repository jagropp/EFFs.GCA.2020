% Calculate the isotopic fractionations for methanogenesis and AOM, the
% different scenarios are presented in Table 10 in the paper.

clear
% Choose temperature in oC
Tc        = 25;
% Generate list of beta values
beta_vals = calc_132_betas(Tc);
% Number of simulations where relevant
sims      = 10; 

% Hydrogenotrophic methanogenesis for three scenarios
[klna_CO2_CH4_i,klna_CO2_CH4_ii,klna_CO2_CH4_iii] = HT_MOG(beta_vals,sims);
% AOM
klna_AOM_CH4_CH2   = AOM(beta_vals);
% Acetoclastic methanogenesis
klna_aceto_CH3_CH4 = aceto_CH3(beta_vals); % From methyl group to CH4
klna_aceto_COO_CO2 = aceto_COO(beta_vals); % From carboxyl group to CO2

function [kln_CO2_CH4_i,kln_CO2_CH4_ii,kln_CO2_CH4_iii] = ...
    HT_MOG(beta_vals,sims)
% Local function to calculate isotopic fractioantion in hydrogenotrophic
% methanogenesis as described in table 7 and in the discussion.

%%(i) WANG model
a13eq(1) = beta_vals.beta13_values(2)./beta_vals.beta13_values(10);
a13eq(2) = beta_vals.beta13_values(10)./beta_vals.beta13_values(11);
a13eq(3) = beta_vals.beta13_values(11)./beta_vals.beta13_values(13);
a13eq(4) = beta_vals.beta13_values(13)./beta_vals.beta13_values(4);

a13kin_f    = ones(1,3).*1.02;
a13kin_f(4) = 1.04;
f_vec       = linspace(0,1,sims);

aCO2_CH4_i = ones(1,length(f_vec));
for sims = 1:length(f_vec)
    for i = fliplr(1:4)
        aCO2_CH4_i(sims) = ...
            (a13eq(i).*aCO2_CH4_i(sims) - a13kin_f(i)).*f_vec(sims) + a13kin_f(i);
    end
end
kln_CO2_CH4_i = round(1000.*log(aCO2_CH4_i),1);

%%(ii) STOLPER model
clear a13eq f_vec a13kin_f

a13eq(1) = beta_vals.beta13_values(2)./beta_vals.beta13_values(14);
a13eq(2) = beta_vals.beta13_values(14)./beta_vals.beta13_values(4);

a13kin_f(1) = 1;
a13kin_f(2) = 1.04;
f_vec       = linspace(0,1,sims);
f_mat(:,1)  = ones(length(f_vec),1);
f_mat(:,2)  = f_vec;

aCO2_CH4_ii = ones(1,length(f_vec));
for sims = 1:length(f_vec)
    for i = [2 1]
        aCO2_CH4_ii(sims) = ...
            (a13eq(i).*aCO2_CH4_ii(sims) - a13kin_f(i)).*f_mat(sims,i) + a13kin_f(i);
    end
end
kln_CO2_CH4_ii = 1000.*log(aCO2_CH4_ii);

%%(iii) CAO 2019 model
a13eq = [1.0139 1.0196 1.0341 1.0020];
a13eq(1) = beta_vals.beta13_values(2)./beta_vals.beta13_values(10);
a13eq(2) = beta_vals.beta13_values(10)./beta_vals.beta13_values(11);
a13eq(3) = beta_vals.beta13_values(11)./beta_vals.beta13_values(13);
a13eq(4) = beta_vals.beta13_values(13)./beta_vals.beta13_values(4);

a13kin_f    = ones(1,3)*1.02;
a13kin_f(4) = 1.04;
f      = [0 0 0 0;
    1 0 0 0;
    1 1 0 0;
    1 1 1 0;
    1 1 1 1];

aCO2_CH4_iii = ones(1,5);
for sims = 1:5
    for i = fliplr(1:4)
        aCO2_CH4_iii(sims) = ...
            (a13eq(i).*aCO2_CH4_iii(sims) - a13kin_f(i)).*f(sims,i) + a13kin_f(i);
    end
end
kln_CO2_CH4_iii = 1000.*log(aCO2_CH4_iii);
end

function klna = AOM(beta_vals)
%% AOM
a13eq(1) = beta_vals.beta13_values(4)./beta_vals.beta13_values(14);
a13eq(2) = beta_vals.beta13_values(14)./beta_vals.beta13_values(13);
a13eq(3) = beta_vals.beta13_values(13)./beta_vals.beta13_values(11);
a13eq(4) = beta_vals.beta13_values(11)./beta_vals.beta13_values(10);
a13eq(5) = beta_vals.beta13_values(10)./beta_vals.beta13_values(9);
a13eq(6) = beta_vals.beta13_values(9)./beta_vals.beta13_values(8);
a13eq(7) = beta_vals.beta13_values(8)./beta_vals.beta13_values(2);

kff_aom(1)   = 1.0387;
% Select the KFFs for reaction 2 to 7. In Table 11 we present results for
% KFF of 5 and 40 permil.
kff_aom(2:7) = ones(1,6).*1.005;
kff_aom(2:7) = ones(1,6).*1.0408;

% Number of increments between 0 and 1 for the reversibility.
p            = 50; 
k            = 1;
f(:,k)       = linspace(1,0,p);

for j = 1:p
    for n = 1:7
        k      = n;
        f      = repmat([1 1 1 1 1 1 1],p,1);
        f(:,k) = linspace(1,0,p);
        anet_aom(j,n) = 1;
        for i = [7 6 5 4 3 2 1]
            anet_aom(j,n) = (a13eq(i).*anet_aom(j,n) - kff_aom(i)).*f(j,i) + kff_aom(i);
        end
    end
end
klna = 1000.*log(anet_aom);

% Plot figure of net isotopic fractionation relative to the reversibility
% of Mcr. 

% figure
% plot(f(:,k),1000.*log(anet_aom))
% set(gca,'xdir','reverse')
% xlabel(sprintf('f(%g)',k))
% ylabel('ln^{13}\alpha CH_4-H_2O')
% legend('Mcr','Mtr','Mer','Mtd','Mch','Ftr','Fmd','Location','northwest')

end

function [knet_aceto] = aceto_CH3(beta_vals)
%% Acetoclastic - methyl
a13eq(1)    = beta_vals.beta13_values(6)./beta_vals.beta13_values(15);
a13eq(2)    = beta_vals.beta13_values(15)./beta_vals.beta13_values(13);
a13eq(3)    = beta_vals.beta13_values(13)./beta_vals.beta13_values(14);
a13eq(4)    = beta_vals.beta13_values(14)./beta_vals.beta13_values(4);
a13kin_f    = ones(1,3)*exp(40/1000);
a13kin_f(4) = 1.03;
% Choose reversibility scenario
f      = [1 1 1 0];

anet_aceto = 1;
for i = [4 3 2 1]
    anet_aceto = (a13eq(i).*anet_aceto - a13kin_f(i)).*f(i) + a13kin_f(i);
end
knet_aceto = 1000.*log(anet_aceto);
end

function [knet_aceto] = aceto_COO(beta_vals)
%% Acetoclastic - CO2
a13eq(1)    = beta_vals.beta13_values(7)./beta_vals.beta13_values(16);
a13eq(2)    = beta_vals.beta13_values(16)./beta_vals.beta13_values(2);
a13kin_f    = ones(1,2)*1.04;
f      = [1 1];

anet_aceto = 1;
for i = [2 1]
    anet_aceto = (a13eq(i).*anet_aceto - a13kin_f(i)).*f(i) + a13kin_f(i);
end
knet_aceto = round(1000.*log(anet_aceto),1);
end