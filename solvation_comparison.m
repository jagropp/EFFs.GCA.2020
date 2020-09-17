%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyse Otake 2008 and Eldridge 2016
%
% 29.12.17
%
% Differences for sulfur isotope fractionation DFT results with an implicit
% and explicit solvation model.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load 'otake2008'
% Remove S8, was only reported for gas phase and therefore not included in
% the comparison
otake08data(4,:) = []; 
spec = {'H_2S','HS^-','S_2','CS_2','SO_2','SO_3','SO_4^2{-}'};
OtabgCoeff33  = otake08data(1:7,1:4);
OtabaqCoeff33 = otake08data(8:14,1:4);
OtabgCoeff34  = otake08data(1:7,5:8);
OtabaqCoeff34 = otake08data(8:14,5:8);
OtabgCoeff36  = otake08data(1:7,9:12);
OtabaqCoeff36 = otake08data(8:14,9:12);

% T = linspace(-70,650,1000) + 273.15; % Original values from the paper
T = linspace(0,100,1000) + 273.15;
for i = 1:7
    Otabg33(i,:) = OtabgCoeff33(i,1).*1e11./(T.^4) + OtabgCoeff33(i,2).*1e8./(T.^3) + ...
        OtabgCoeff33(i,3).*1e6./(T.^2) + OtabgCoeff33(i,4).*1e3./T;
    Otabaq33(i,:) = OtabaqCoeff33(i,1).*1e11./(T.^4) + OtabaqCoeff33(i,2).*1e8./(T.^3) + ...
        OtabaqCoeff33(i,3).*1e6./(T.^2) + OtabaqCoeff33(i,4).*1e3./T;
    Otabg34(i,:) = OtabgCoeff34(i,1).*1e11./(T.^4) + OtabgCoeff34(i,2).*1e8./(T.^3) + ...
        OtabgCoeff34(i,3).*1e6./(T.^2) + OtabgCoeff34(i,4).*1e3./T;
    Otabaq34(i,:) = OtabaqCoeff34(i,1).*1e11./(T.^4) + OtabaqCoeff34(i,2).*1e8./(T.^3) + ...
        OtabaqCoeff34(i,3).*1e6./(T.^2) + OtabaqCoeff34(i,4).*1e3./T;
    Otabg36(i,:) = OtabgCoeff36(i,1).*1e11./(T.^4) + OtabgCoeff36(i,2).*1e8./(T.^3) + ...
        OtabgCoeff36(i,3).*1e6./(T.^2) + OtabgCoeff36(i,4).*1e3./T;
    Otabaq36(i,:) = OtabaqCoeff36(i,1).*1e11./(T.^4) + OtabaqCoeff36(i,2).*1e8./(T.^3) + ...
        OtabaqCoeff36(i,3).*1e6./(T.^2) + OtabaqCoeff36(i,4).*1e3./T;
end

plot(1e6./(T.^2),Otabg33)
legend(spec)
Otad33 = abs(Otabg33 - Otabaq33);
Otad34 = abs(Otabg34 - Otabaq34);
Otad36 = abs(Otabg36 - Otabaq36);

figure
subplot(1,3,1)
% line([1e9./((273.15+60).^3) 1e9./((273.15+60).^3)],[0 3])
hold on
plot(1e9./(T.^3),Otad33)
legend(spec)

subplot(1,3,2)
% line([1e9./((273.15+60).^3) 1e9./((273.15+60).^3)],[0 3])
hold on
plot(1e9./(T.^3),(Otad34./1000)+1)

subplot(1,3,3)
% line([1e9./((273.15+60).^3) 1e9./((273.15+60).^3)],[0 3])
hold on
plot(1e9./(T.^3),Otad36)

%%

load 'eldridge2016'
EldbCoeff34  = eldridge2016data(1:28,2:6);
EldbCoeff34(:,1) = EldbCoeff34(:,1)./1e-5;
EldbCoeff34(:,2) = EldbCoeff34(:,2)./1e-4;
EldbCoeff34(:,3) = EldbCoeff34(:,3)./1e-2;
EldbCoeff34(:,4) = EldbCoeff34(:,4)./1e2;

for i = 1:28
    Eldb34(i,:) = EldbCoeff34(i,1)./(T.^4) + EldbCoeff34(i,2)./(T.^3) + EldbCoeff34(i,3)./(T.^2) + ...
                     EldbCoeff34(i,4)./T + EldbCoeff34(i,5);
end

% eld = [1 4 8 9 15 17 18 19]; % both 4 and 5 are SO3(2-)
eld = [1 5 9 10 16 18 19 20];
ota = [14 13 12 5 8 9 1 10];
Otab34 = [Otabg34; Otabaq34];
spec1 = {'SO_4^{2-}','SO_3','SO_{2(aq)}','SO_{2(g)}','H_2S_{(aq)}','HS^-','H_2S_{(g)}','S^{2-}'};
dif = (1000.*(Eldb34(eld,:)-1) - (Otab34(ota,:)));

%%
figure

subplot(1,2,1)
for i = 1:length(eld)
    plot(1000./T,(Otab34(ota(i),:)./1000)+1)
    text(1000./T(1)+0.015,(Otab34(ota(i),1)./1000)+1,spec1{i},'FontName','Myriad pro')
    hold on
end
% title('Otake (2008) - Implicit solvation model')
ylabel('^{34}\beta')
xlabel('1000/T [K]')
grid on
grid minor

subplot(1,2,2)
for i = 1:length(eld)
    plot(1000./T,Eldb34(eld(i),:))
    text(1000./T(1)+0.015,Eldb34(eld(i),1),spec1{i},'FontName','Myriad pro')
    hold on
end
% title('Eldbridge (2017) - Explicit solvation model')
ylabel('^{34}\beta')
xlabel('1000/T [K]')
grid on
grid minor

%%
figure
for i = 1:length(eld)
    plot(1000./T,dif(i,:))
    text(1000./T(1)+0.015,dif(i),spec1{i},'FontName','Myriad pro')
    hold on
end
ylabel('\Delta^{34}\beta_{Eldbridge-Otake}')
xlabel('1000/T [K]')
grid on
grid minor

















