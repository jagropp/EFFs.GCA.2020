% 
% maxT = 400;
% T = linspace(0,maxT,maxT);
% T = T + 273.15;

% % 1.  H2O(l) <-> H2O(g) - Horita 1994 T = 0-374 oc (experimental)
% % 2.  H2O(g) <-> H2(g) - Bottinga 1969 T = 0-600 oc
% % 3.  H2O(g) <-> H2(g) - Suess 1949 T = 0-1000 oc
% 4.  H2O(l) <-> H2(g) - Rolston 1976 T = 0-200 oc
% % 5.  CH4(g) <-> H2(g) - Bottinga 1969 T = 0-700 oc
% % 6.  H2O(g) <-> CH4(g) - Bottinga 1969 T = 10-250 oc
% % 7.  CH4(g) <-> H2(g) - Horibe & Craig 1995 T = 200-500 oc
% % 8.  H2O(l) <-> H2O(g) - Majoube 1971 T = 0-100 oc (theoretical)
% % 9.  H2O(l) <-> H2O(g) - Merlivat 1967 T = -15 - 0 oc (theoretical)
% % 10. H2O(l) <-> H2O(g) - Kakiuchi 1979 T = 10-40 oc (experimental)
% 11. CH4(g) <-> H2O(l) - calculation -> Horibe - Rolston (7-4)
% 12. CH4(g) <-> H2O(l) - calculation -> -(Horita + Bottinga) -(1+6)
% 13. H2O(g) <-> CH4(g) - calculation -> 2 + (-5)
% 14. H2O(g) <-> CH4(g) - calculation -> 2 + (-7)
% 15. H2O(g) <-> CH4(g) - calculation -> 3 + (-5)
% 16. H2O(g) <-> CH4(g) - calculation -> 3 + (-7)
% 17. H2O(g) <-> CH4(g) - Mark Iron, 17.12.18 M06L TZVP W06
% 18. H2O(l) <-> H2(g) - Horibe & Craig 1995 T = 0 - 370 oc (eqn. 8)
% 19. H2O(l) <-> CH4(g) - Horibe - Horibe

A = [  -0.353 0 0 0 0 0 0 0 0 0];
B = [  42.170 0 0 0 0 0 0 0 0 0];
C = [-309.400 0 0 0 0 0 0 0 0 0];
D = [  963.700   13.000    0.000   27.870   25.000 -7.690 0  24.844   15.013    2.408]; 
E = [-1399.000  389.610  467.600  368.900  346.000  6.100 0 -76.248    0.000   64.550];
F = [  766.200 -204.340 -304.000 -214.300 -223.000 88.400 0  56.612 -100.000 -168.000];

% Original values are in comment
Tval = [0 400; % 0-374
        0 400; % 0-600
        0 400; % 0-1000
        0 200;
        0 400; % 0-700
        10 250; % 10-250
        0 500; % 200-500
        0 100;
        -15 0;
        10 40;
        0 400;
        0 400;
        0 400;
        0 400;
        0 400;
        0 400;
        0 700;
        0 400;
        0 400];
Tval = Tval + 273.15;

vals = 400;
siz = size(Tval);
for i = 1:siz(1)
    T(i,:) = linspace(Tval(i,1),Tval(i,2),vals);
end

eps = zeros(length(A),length(T));
for i = 1:length(A)
eps(i,:) = A(i).*(1e18./T(i,:).^6) + B(i).*(1e12./T(i,:).^4) + C(i).*(1e9./T(i,:).^3) + ...
            D(i).*(1e6./T(i,:).^2) + E(i).*(1e3./T(i,:)) + F(i); % 1000 ln alpha
end

eps(1,:) = 1158.8.*(T(1,:).^3./1e9) - 1620.*(T(1,:).^2./1e6) + 794.84.*(T(1,:)./1e3) -...
            161.04 + 2.992.*(1e9./T(1,:).^3);
% eps(1,373:maxT) = NaN;

% A correction to get H2O(l) <-> H2(g), as was prented in Wang 2015 (sup).
eps(2,:) = eps(2,:) + eps(1,:);
eps(3,:) = eps(3,:) + eps(1,:);

% Values from Horibe and Craig, GCA, 1995. The equation here was presented
% in the original paper and there it refers to alpha.
eps(7,:) = 1000.*log(0.8994 + 183.540e3./(T(7,:).^2));
% Calculated value of CH4(g) <-> H2O(l) - same as in Wang
eps(11,:) = eps(7,:) - eps(4,:);
% Another way to calculate CH4(g) <-> H2O(l) - 
eps(12,:) = -((eps(1,:) + eps(6,:)));
% H2O(g) <-> CH4(g) - calculation
eps(13,:) = eps(2,:) - eps(5,:);
eps(14,:) = eps(2,:) - eps(7,:);
eps(15,:) = eps(3,:) - eps(5,:);
eps(16,:) = eps(3,:) - eps(7,:);
% Mark Iron
% eps(17,:) = 9e-13.*T(17,:).^5 - 4e-9.*T(17,:).^4 + 5e-6.*T(17,:).^3 ... 
%                 - 0.0044.*T(17,:).^2 + 1.7995.*T(17,:) + 32.946;                
eps(17,:) = 7e-13.*T(17,:).^5 - 3e-09.*T(17,:).^4 + 5e-06.*T(17,:).^3 - 0.0037.*T(17,:).^2 + 1.4634.*T(17,:) - 85.238;
eps(18,:) = 1000.*log(1.0473 + 201036./T(18,:).^2 + 2.06e9./T(18,:).^4 + 0.18e15./T(18,:).^6);
eps(19,:) = eps(18,:) - eps(7,:);
%%

% Fig. S3 from Wang et al. 2015.
figure(1);
set(1,'Units','centimeters','Position',[10 5 12 10])
pos = [2 3 4 7 18 19];
for i = 1:length(pos)
    plot(T(pos(i),:)-273.15,exp(eps(pos(i),:)./1000))
%     plot(T(pos(i),:)-273.15,eps(pos(i),:))
    hold on
end
xlim([0 400])
box off
xlabel('T [^oC]')
% xlabel('10^6/T [K]')
ylabel('\alpha')
l = legend(...
       'H$_2$O(l) - H$_2$(g) Bottinga 1969 $\times$ Horita 1994',...
       'H$_2$O(l) - H$_2$(g) Suess 1949 $\times$ Horita 1994',...
       'H$_2$O(l) - H$_2$(g) Rolston 1976',...
       'CH$_4$(g) - H$_2$(g) Horibe \& Craig 1995',...
       'H$_2$O(l) - H$_2$(g) Horibe \& Craig 1995',...
       'H$_2$O(l) - CH$_4$(g) Horibe \& Craig 1995',...
       'Location','best');
set(l,'FontSize',10)
set(l,'interpreter','latex')
legend('boxoff')

% H2O(l) <-> H2O(g) compilation
figure(2);
set(2,'Units','centimeters','Position',[10 5 12 10])
pos = [1 8 9 10];
for i = 1:length(pos)
    plot(T(pos(i),:)-273.15,exp(eps(pos(i),:)./1000))
    hold on
end
xlim([-50 400])
box off
xlabel('T [^oC]')
ylabel('\alpha')
l = legend('Horita 1994; T = 0-374$^o$C',...
           'Majoube 1971; T = 0-100$^o$C',...
           'Merlivat 1967; T = -15-0$^o$C',...
           'Kakiuchi 1979; T = 10-40$^o$C',...
           'Location','best');
set(l,'FontSize',10)
set(l,'interpreter','latex')
legend('boxoff')
%% H2O(g) <-> CH4(g) compilation
figure(3);
set(3,'Units','centimeters','Position',[10 5 12 10])
pos = [6 13 14 15 16 17];
line = {'-';'-k';'--k';'-.k';':k';'--r'};
for i = 1:length(pos)
    plot(T(pos(i),:)-273.15,exp(eps(pos(i),:)./1000),line{i})
    hold on
end
% xlim([0 400])
% ylim([0.5 4])
box off
xlabel('T [^oC]')
ylabel('\alpha')
l = legend('H$_2$O(g) - CH$_4$(g) Bottinga 1969',... 
           'H$_2$O(g) - CH$_4$(g) B69 $\times$ B69',... 
           'H$_2$O(g) - CH$_4$(g) B69 $\times$ H95',... 
           'H$_2$O(g) - CH$_4$(g) S49 $\times$ B69',... 
           'H$_2$O(g) - CH$_4$(g) S49 $\times$ H95',... 
           'H$_2$O(g) - CH$_4$(g) M. Iron HCTH',... 
           'Location','best');
set(l,'FontSize',10)
set(l,'interpreter','latex')
legend('boxoff')
xlim([0 400])

%% H2O(g) <-> CH4(g) compilation
figure(4);
set(4,'Units','centimeters','Position',[10 5 12 10])
pos = [2 3 5 7];
% line = {'-';'-';'-';'-k';'--k';'-.k';':k';'--r'};
for i = 1:length(pos)-1
%     plot(T(pos(i),:)-273.15,exp(eps(pos(i),:)./1000),line{i})
    plot(T(pos(i),:)-273.15,exp(eps(pos(i),:)./1000))    
    hold on
end
ind = find(T(7,:)>=200+273.15);
plot(T(7,ind)-273.15,exp(eps(7,ind)./1000))    
xlim([0 500])
% ylim([0.5 4])
box off
xlabel('T [^oC]')
ylabel('\alpha')
l = legend('H$_2$O(g) - H$_2$(g) Bottinga 1969',...
           'H$_2$O(g) - H$_2$(g) Suess 1949',...
           'CH$_4$(g) - H$_2$(g) Bottinga 1969',...
           'CH$_4$(g) - H$_2$(g) Horibe \& Craig 1995',...
           'Location','best');
set(l,'FontSize',10)
set(l,'interpreter','latex')
legend('boxoff')
