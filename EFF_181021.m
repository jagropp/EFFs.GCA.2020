% 
maxT = 701;
Tc = linspace(0,maxT-1,maxT);
Tk = Tc + 273.15;
[b1,b2] = xlsread('C:\Jonathan\Dropbox (Weizmann Institute)\Jonathan_(with_Itay)\Paper I EFFs\alpha-beta.xlsx',6);
b1(1,:) = [];
Tc_DFT = b1(:,1);
Tk_DFT = Tc_DFT + 273.15;

%% Hydrogen

% 1. H2O(l) <-> H2O(g) - Horita 1994 T = 0-374 oc (Experimental)
eps(1,:) = 1158.8.*(Tk.^3./1e9) - 1620.1.*(Tk.^2./1e6) + 794.84.*(Tk./1e3) - 161.04 + 2.9992.*(1e9./Tk.^3);
eps(1,375:end) = nan;
eps_h2o_lg_DFT = 1158.8.*(Tk_DFT.^3./1e9) - 1620.1.*(Tk_DFT.^2./1e6) + 794.84.*(Tk_DFT./1e3) - 161.04 + 2.9992.*(1e9./Tk_DFT.^3);

% 2. H2O(l) <-> H2O(g) - Majoube 1971 T = 0-100 oc (Theoretical)
eps(2,:) = 24.844.*(1e6./Tk.^2) - 76.248.*(1e3./Tk) + 52.612;
eps(2,101:end) = nan;

% 3. H2O(l) <-> H2O(g) - Kakiuchi 1979 T = 10-40 oc (Experimental)
eps(3,:) = 2.408.*(1e6./Tk.^2) + 64.55.*(1e3./Tk) - 168;
eps(3,[1:9,41:end]) = nan;

% 4. H2O(g) <-> H2(g) - Bottinga 1969 T = 0-600 oc (Theoretical)
eps(4,:) = 13.*(1e6./Tk.^2) + 389.61.*(1e3./Tk) - 204.34;
eps(4,601:end) = nan;

% 5. H2O(g) <-> H2(g) - Suess 1949 T = 0-1000 oc (Theoretical)
eps(5,:) = 467.6.*(1e3./Tk) - 304;

% 6. H2O(g) <-> H2(g) - Cerrai 1954 T = 50-742 oc (Experimental)
eps(6,:) = 1000.*log(1.045 + 211286./(Tk.^2)); % Eq. from Horibe 1995
eps(6,1:49) = nan;

% 7. H2O(l) <-> H2(g) - Rolston 1976 T = 0-200 oc (Experimental)
eps(7,:) = 27.87.*(1e6./Tk.^2) + 368.9.*(1e3./Tk) - 214.3;
eps(7,201:end) = nan;

% 8. H2O(l) <-> H2(g) - Horibe & Craig 1995 T = 0 - 370 oc (Theoretical)
eps(8,:) = 1000.*log(1.0473 + 201036./Tk.^2 + 2.06.*1e9./Tk.^4 + 0.18.*1e15./Tk.^6); % Eqn. 8 in the paper

% 9. CH4(g) <-> H2(g) - Bottinga 1969 T = 0-700 oc (Theoretical)
eps(9,:) = 25.*(1e6./Tk.^2) + 346.*(1e3./Tk) -223;

% 10. H2O(g) <-> CH4(g) - Bottinga 1969 T = 10-250 oc (Theoretical)
eps(10,:) = -7.69.*(1e6./Tk.^2) + 6.1.*(1e3./Tk) + 88.4;
eps(10,[1:9,251:end]) = nan;
% 11.  CH4(g) <-> H2(g) - Horibe & Craig 1995 T = 200-500 oc (Experimental)
eps(11,:) = 1000.*log(0.8994 + 183540./Tk.^2);
eps(11,[1:199,501:end]) = nan;

% eps(i,:) = A.*(1e18./Tk.^6) + B.*(1e12./Tk.^4) + C.*(1e9./Tk.^3) + ...
%            D.*(1e6./Tk.^2) + E.*(1e3./Tk) + F;

%% Combinations of fractionation factors

% H2O(l)-CH4(g) based on H2O(l)-H2O(g) from Horita 
eps(12,:) = eps(1,:) + eps(10,:); % Bottinga 1
eps(13,:) = eps(1,:) + eps(4,:) - eps(9,:); % Bottinga 2
eps(14,:) = eps(1,:) + eps(5,:) - eps(11,:); % Suess 49 + Horibe 95
eps(15,:) = eps(1,:) + eps(6,:) - eps(11,:); % Cerrai 49 + Horibe 95
% eps(16,:) = eps(7,:) + eps(5,:) - eps(11,:); % Suess 49 + Horibe 95
% eps(17,:) = eps(7,:) + eps(6,:) - eps(11,:); % Cerrai 49 + Horibe 95

%% H2O(l) <-> CH4(g) compilation
% based on H2O(l <-> g) of Horita
figure(1);
set(1,'Units','centimeters','Position',[10 5 10 8])
pos = 12:15;
line = {'-k';':k';'--k';'-.k'};
for i = 1:length(pos)
    plot(Tc,eps(pos(i),:),line{i})
    hold on
end
plot(b1(:,1),eps_h2o_lg_DFT+b1(:,29))
xlim([0 374])
box off
xlabel('T [\circC]')
ylabel('1000\timesln^2\alpha')
l = legend('B69','B69\timesB69','S49\timesH95','C49\timesH95','This work','Location','best');
set(l,'FontSize',10)
legend('boxoff')

%% H2O(g) - H2(g)