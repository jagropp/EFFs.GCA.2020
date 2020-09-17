% Plot predicted CNO isotope equilibrium fractionation factors (EFFs)
% compated to experimental results.

[dat,lab] = xlsread('181122_CNO_EIF.xlsx');
for i = 1:8; defs{i} = lab{1,i+3}; end

expdat = 1000.*log(dat(1:37,4));
dftdat = 1000.*log(dat(1:37,5:12));
expdat([21 22 26 32 33 34])   = nan; 
dftdat([21 22 26 32 33 34],:) = nan;
ccode = [255,247,251
    236,231,242
    208,209,230
    166,189,219
    116,169,207
    54,144,192
    5,112,176
    3,78,123]./255;
%%
figure(1);
set(1,'Units','centimeters','Position',[10 10 9 8])
ha1 = line([-20 40],[-20 40]);
ha1.Color = [1 1 1].*0.5;
ha1.LineWidth = 0.5;
hold on
for ifig = 1:8
    ha2 = scatter(expdat,dftdat(:,ifig));
    ha2.MarkerEdgeColor = 'k';
    ha2.MarkerFaceColor = ccode(ifig,:);
    ha2.LineWidth = 0.5;
    hha(ifig) = ha2;
    hold on
end
l = legend([hha(1) hha(2) hha(3) hha(4)  hha(5) hha(6) hha(7) hha(8)],defs);
l.Box = 'off';
l.Location = 'northwest';
xlim([-20 40]);
ylim([-20 40]);
set(gca,'YTick',[-10 0 10 20 30])
xlabel('1000ln\alpha (Experimental)')
ylabel('1000ln\alpha (Prediction)')

%% 
ccode2 = [31,120,180
178,223,138
51,160,44]./255;
% ccode2 = [0 0.447 0.741; 0.850 0.325 0.098; 0.929 0.694 0.125];
% 95% predicition interval
rang = [1 20; 21 27; 28 37];
x = expdat;
y = dftdat(:,1);
x(isnan(x)) = [];
y(isnan(y)) = [];
[fitresult,SSE,R2,stdError] = fit(x,y,'poly1');
ci = predint(fitresult,sort(x),0.95,'Functional','on'); % Simultaneous Functional Bounds

xC = expdat(1:20);  yC = dftdat(1:20);
xN = expdat(21:27); yN = dftdat(21:27);
xO = expdat(28:37); yO = dftdat(28:37);
xC(isnan(xC)) = []; yC(isnan(yC)) = [];
xN(isnan(xN)) = []; yN(isnan(yN)) = [];
xO(isnan(xO)) = []; yO(isnan(yO)) = [];
[fitresultC,SSEC,~,stdErrorC] = fit(xC,yC','poly1');
[fitresultN,SSEN,~,stdErrorN] = fit(xN,yN','poly1');
[fitresultO,SSEO,~,stdErrorO] = fit(xO,yO','poly1');

figure(2);
set(2,'Units','centimeters','Position',[10 10 9 8])
hb1 = line([-20 40],[-20 40]);
hb1.Color = [1 1 1].*0;
hb1.LineWidth = 1;
hold on
hb2 = plot(sort(x),ci(:,1),':',sort(x),ci(:,2),':','Color',[1 1 1].*0.5);
hold on
hb3 = plot(sort(x),fitresult.p1.*sort(x)+fitresult.p2,'Color',[1 1 1].*0.5);
hold on
for ifig = 1:3
    hb4 = scatter(expdat(rang(ifig,1):rang(ifig,2)),dftdat((rang(ifig,1):rang(ifig,2)),1));
    hb4.MarkerEdgeColor = 'k';
    hb4.MarkerFaceColor = ccode2(ifig,:);
    hb4.LineWidth = 0.5;
    hhb(ifig) = hb4;
end

if fitresult.p2 > 0
    tab = text(12,-5,sprintf('y = %.3fx + %.3f',fitresult.p1,fitresult.p2));
elseif fitresult.p2 <= 0
    tab = text(12,-5,sprintf('y = %.3fx - %.3f',fitresult.p1,-fitresult.p2));
end
tab.FontName = 'Acumin pro';
tR2 = text(12,-9,-160,sprintf('R^2 = %.3f',SSE.rsquare));
tR2.FontName = 'Acumin pro';

l = legend([hhb(1) hhb(2) hhb(3)],{'C','N','O'});
l.Box = 'off';
l.Location = 'northwest';
xlim([-15 35]);
ylim([-15 35]);
set(gca,'YTick',[-10 0 10 20 30])
xlabel('1000ln\alpha (Experimental)')
ylabel('1000ln\alpha (Prediction)')

%%
pathstr = 'C:\Jonathan\Dropbox (Weizmann Institute)\Group meetings\181128_Jonathan\';
print('-f1',sprintf('%sCNO',pathstr),'-dpng','-r600');
print('-f1',sprintf('%sCNO',pathstr),'-dsvg');
print('-f2',sprintf('%sCNO_M06L',pathstr),'-dpng','-r600');
print('-f2',sprintf('%sCNO_M06L',pathstr),'-dsvg');

%%
% Export sheet names in 181121 file
filepath = 'C:\Jonathan\Dropbox (Weizmann Institute)\Collaborations\Mark Iron\181121 Equilibrium Isotopic Fractionation - Liu';
xlRange = 'D26:D38';
[status,sheets] = xlsfinfo(filepath);
l = 48;
dati = zeros(12,l);
for i = 1:l
    dati(:,i) = xlsread(filepath,i,xlRange);
end
datidft = 1000.*(dati - 1);
datiexp(:,1:28)  = repmat(xlsread(filepath,1,'E26:E38'),1,28);
datiexp(:,29:38) = repmat(xlsread(filepath,29,'E26:E38'),1,10);
datiexp(:,39:48) = repmat(xlsread(filepath,39,'E26:E38'),1,10);
datiexp(11,:) = [];
datidft(11,:) = [];
%%
mae = abs(datidft-datiexp);
mae(end+1,:) = nanmean(mae,1);
[~,mae25pos] = sort(mae(end,1:28));
[~,mae50pos] = sort(mae(end,29:38));
[~,mae70pos] = sort(mae(end,39:48));

clear xi yi
lin = [2 3 4 6 8 10];
cyc = [1 5 7 9 11];
all = 1:11;
lincycrange = cyc;
filenameext = '_cyc';
xi = datiexp(lincycrange,[1 29 39]);
yi = datidft(lincycrange,[mae25pos(1) mae50pos(1)+28 mae70pos(1)+38]);
% xi(isnan(xi)) = [];
% yi(isnan(yi)) = [];
xi = reshape(xi,[length(xi)*3,1]);
yi = reshape(yi,[length(yi)*3,1]);
[xii,posxi] = sort(xi);
yii = yi(posxi);
[fitresult,SSE,R2,stdError] = fit(xii,yii,'poly1');
ciHD = predint(fitresult,xii,0.95,'Functional','on'); % Simultaneous Functional Bounds

figure(3);
set(3,'Units','centimeters','Position',[10 10 9 8])
hc1 = line([-300 100],[-300 100],'Color','k');
hold on
hc2 = plot(xii,ciHD(:,1),':',xii,ciHD(:,2),':','Color',[1 1 1].*0.5);
hold on
hc3 = plot(xii,fitresult.p1.*xii+fitresult.p2,'Color',[1 1 1].*0.5);
hold on
hc4 = plot(datiexp(lincycrange,1),datidft(lincycrange,mae25pos(1)),'o');
hc4.MarkerFaceColor = [0.000 0.447 0.741];
hc4.MarkerEdgeColor = 'k';
hc4.LineWidth = 0.5;
hold on
hc5 = plot(datiexp(lincycrange,29),datidft(lincycrange,mae50pos(1)+28),'o');
hc5.MarkerFaceColor = [0.850 0.325 0.098];
hc5.MarkerEdgeColor = 'k';
hc5.LineWidth = 0.5;
hold on
hc6 = plot(datiexp(lincycrange,39),datidft(lincycrange,mae70pos(1)+38),'o');
hc6.MarkerFaceColor = [0.929 0.694 0.125];
hc6.MarkerEdgeColor = 'k';
hc6.LineWidth = 0.5;

if fitresult.p2 > 0
    tab = text(-70,-140,sprintf('y = %.3fx + %.3f',fitresult.p1,fitresult.p2));
elseif fitresult.p2 <= 0
    tab = text(-70,-140,sprintf('y = %.3fx - %.3f',fitresult.p1,-fitresult.p2));
end
tab.FontName = 'Acumin pro';
tR2 = text(-70,-160,sprintf('R^2 = %.3f',SSE.rsquare));
tR2.FontName = 'Acumin pro';

l = legend([hc4 hc5 hc6],{'25\circC','50\circC','70\circC'});
l.Location = 'northwest';
l.Box = 'off';
xlim([-200 50])
ylim([-200 50])
xlabel('1000ln^2\alpha (Experimental)')
ylabel('1000ln^2\alpha (Prediction)')

print('-f3',sprintf('%sHD-Ketones%s',pathstr,filenameext),'-dpng','-r600');
print('-f3',sprintf('%sHD-Ketones%s',pathstr,filenameext),'-dsvg');




