clear
fldr = pwd;
% addpath('C:\Jonathan\Dropbox (Weizmann Institute)\Collaborations\Mark Iron')
cd('/Users/jonathag/Dropbox (Weizmann Institute)/Collaborations/Mark Iron')
[dat,lab] = xlsread('190430_exp_calc_data',1);
cd(fldr)
atom   = char(lab{2:end,1});
oliers = char(lab{2:end,9});
keto   = char(lab{2:end,10});

c = [178,24,43
    214,96,77
    244,165,130
    253,219,199
    247,247,247
    209,229,240
    146,197,222
    67,147,195
    33,102,172]./255;

ind = 35;
xdat = dat(1:ind,2);
clear ydat
ydat = dat(1:ind,6); % M06-L D3

xdatHlin = dat(keto=='L',2);
ydatHlin = dat(keto=='L',6); % H M06-L D3 TZVP
xdatHcyc = dat(keto=='C',2);
ydatHcyc = dat(keto=='C',6); % H M06-L D3 TZVP
xdatH    = dat(keto=='C' | keto=='L',2);
ydatH    = dat(keto=='C' | keto=='L',6);
%%
% Remove ouliers from calculation of CI and PI, by assigning a threshold
% (trsh) in units of permil.
err = 1000.*abs(xdat-ydat);
trsh = 10;
% trsh = 30;
xdat_i = xdat;
ydat_i = ydat;
atom_i = atom;
xdat_i(err>trsh) = [];
ydat_i(err>trsh) = [];
atom_i(err>trsh) = [];

cCNO                = zeros(length(xdat_i),3);
cCNO(atom_i == 'C',:) = repmat(c(1,:),length(xdat_i(atom_i == 'C')),1);
cCNO(atom_i == 'N',:) = repmat(c(5,:),length(xdat_i(atom_i == 'N')),1);
cCNO(atom_i == 'O',:) = repmat(c(9,:),length(xdat_i(atom_i == 'O')),1);

lim = [min([min(xdat_i) min(ydat_i)])-0.005 max([max(xdat_i) max(ydat_i)])+0.005];

xfit    = linspace(lim(1),lim(2),1000);
m0      = [1 0];
[s1,s2] = fit(xdat_i,ydat_i,'poly1');
m       = [s1.p1 s1.p2];
yfit    = m(1).*xfit + m(2);
rsq     = s2.rsquare;
% PREDICITION INTERVAL - HERE IT IS 95%!!! i.e., 2 sigma
pilevel = 0.95;
PI95    = predint(s1,xfit,pilevel,'observation','off');

[sa1,sa2] = polyfit(xdat_i,ydat_i,1);
[~,dyCI] = polyconf(sa1,xfit,sa2,'predopt','curve','Alpha',1-pilevel);

marksize = 18;
col = 'k';
mark = 'o';
defFont = {'Royal Society of Chemistry'};
% def
tYpos = lim(2)-0.003;

h = figure;
set(h,'Units','Centimeters','Position',[2 2 8.3 7])

% Plot linear regression
plot(xfit,yfit,col,'LineWidth',0.5); hold on
% Plot predicition interval
plot(xfit,PI95,':','Color',col,'LineWidth',0.5); hold on
% Plot confidence interval
plot(xfit,yfit+dyCI,'--','Color','k','LineWidth',0.5); hold on
plot(xfit,yfit-dyCI,'--','Color','k','LineWidth',0.5); hold on

set(gca,'FontSize',8)

for j = 1:length(xdat_i)
    h1 = scatter(xdat_i(j),ydat_i(j),marksize,'filled');
    h1.MarkerFaceColor = cCNO(j,:);
    h1.Marker = mark;
    h1.MarkerEdgeColor = col;
end
hold on
indErr = find(err>trsh);
for j = 1:length(indErr)
    h2 = scatter(xdat(indErr(j)),ydat(indErr(j)));
    h2.Marker = '+';
    h2.MarkerEdgeColor = 'k';
end
if m(2) > 0
    t = text(lim(1)+0.003,tYpos,sprintf('y = %.3fx + %.3f, R^2 = %.2f',m(1),m(2),rsq));
elseif m(2) < 0
    t = text(lim(1)+0.003,tYpos,sprintf('y = %.3fx - %.3f, R^2 = %.2f',m(1),abs(m(2)),rsq));
end
% t = text(lim(1)+0.003,tYpos,'y = 0.9975x, R^2 = 0.77');
t.FontName = defFont{1};
t.FontSize = 7.5;
t.Color = col;
hold on

f = get(gca,'Children');
f = flipud(f);
ind(1) = find(cCNO(:,1)==c(1,1),1,'first')+5;
ind(2) = find(cCNO(:,1)==c(5,1),1,'first')+5;
ind(3) = find(cCNO(:,1)==c(9,1),1,'first')+5;
l = legend([f(ind(1)) f(ind(2)) f(ind(3))],{'C','N','O'});
l.Location = 'southeast';
l.Box = 'on';
% l.EdgeColor = 'w';

axTicks  = [0.99,1.00,1.01,1.02,1.03];
axLabels = {'0.99','1.00','1.01','1.02','1.03'};
set(gca,'XTickLabel',axLabels,'YTickLabel',axLabels,...
    'XTick',axTicks,'YTick',axTicks)
xl = xlabel('EIF (Experimental)','FontSize',10);
xl.FontName = defFont{1};
yl = ylabel('EIF (Calculated)','FontSize',10);
yl.FontName = defFont{1};

grid on
xlim(lim)
ylim(lim)

%%
limH = [min([min([xdatHlin;xdatHcyc]) min([ydatHlin;ydatHcyc])])-0.03...
        max([max([xdatHlin;xdatHcyc]) max([ydatHlin;ydatHcyc])])+0.03];
xfitH = linspace(limH(1),limH(2),1000);
% 
% % cond = 'L';
% % xdatH(keto==cond) = [];
% % ydatH(keto==cond) = [];
% % keto(keto==cond)  = [];
% 
[s1lin,s2lin] = fit(xdatHlin,ydatHlin,'poly1');
mlin       = [s1lin.p1 s1lin.p2];
yfitHlin   = mlin(1).*xfitH + mlin(2);
rsqlin     = s2lin.rsquare;
PI95lin    = predint(s1lin,xfitH,0.95,'observation','off');
[sa1lin,sa2lin] = polyfit(xdatHlin,ydatHlin,1);
[~,dyCIlin] = polyconf(sa1lin,xfitH,sa2lin,'predopt','curve','Alpha',1-pilevel);

[s1cyc,s2cyc] = fit(xdatHcyc,ydatHcyc,'poly1');
mcyc       = [s1cyc.p1 s1cyc.p2];
yfitHcyc   = mcyc(1).*xfitH + mcyc(2);
rsqcyc     = s2cyc.rsquare;
PI95cyc    = predint(s1cyc,xfitH,0.95,'observation','off');
[sa1cyc,sa2cyc] = polyfit(xdatHcyc,ydatHcyc,1);
[~,dyCIcyc] = polyconf(sa1cyc,xfitH,sa2cyc,'predopt','curve','Alpha',1-pilevel);
% 
% marksize = 18;
% defFont = {'Royal Society of Chemistry'};
% tYpos = limH(2)-0.02;
% 
% linCol = [166,217,106]./255;
% cycCol = [253,174,97]./255;
% 
% h = figure;
% set(h,'Units','Centimeters','Position',[2 2 8.3 7])
% 
% plot(xfitH,yfitHlin,'Color',linCol,'LineWidth',1.5); hold on
% plot(xfitH,PI95lin,':','Color',linCol,'LineWidth',1.5); hold on
% 
% plot(xfitH,yfitHlin+dyCIlin,'-','Color',c(1,:),'LineWidth',0.5); hold on
% plot(xfitH,yfitHlin-dyCIlin,'-','Color',c(1,:),'LineWidth',0.5); hold on
% 
% plot(xfitH,yfitHcyc,'Color',cycCol,'LineWidth',1.5); hold on
% plot(xfitH,PI95cyc,':','Color',cycCol,'LineWidth',1.5); hold on
% 
% set(gca,'FontSize',8)
% 
% for j = 1:length(xdatHlin)
%     h1 = scatter(xdatHlin(j),ydatHlin(j),marksize,'filled');
%     h1.MarkerFaceColor = linCol;
%     h1.MarkerEdgeColor = 'k';
% end
% 
% for j = 1:length(xdatHcyc)
%     h1 = scatter(xdatHcyc(j),ydatHcyc(j),marksize,'filled');
%     h1.MarkerFaceColor = cycCol;
%     h1.MarkerEdgeColor = 'k';
% end
% 
% % if m(2) > 0
% %     t = text(limH(1)+0.02,tYpos,sprintf('y = %.3fx + %.3f, R^2 = %.2f',m(1),m(2),rsq));
% % elseif m(2) < 0
% %     t = text(limH(1)+0.02,tYpos,sprintf('y = %.3fx - %.3f, R^2 = %.2f',m(1),abs(m(2)),rsq));
% % end
% % t.FontName = defFont{1};
% % t.FontSize = 7.5;
% % t.Color = 'k';
% hold on
% 
% f = get(gca,'Children');
% f = flipud(f);
% % l = legend([f(19) f(4)],{'Linear','Cyclic'});
% % l.Location = 'southeast';
% % l.Box = 'on';
% % l.EdgeColor = 'w';
% 
% axLabels = {'0.85','0.90','0.95','1.00'};
% set(gca,'XTickLabel',axLabels,'YTickLabel',axLabels)
% xl = xlabel('\alpha (Experimental)','FontSize',10);
% xl.FontName = defFont{1};
% yl = ylabel('\alpha (Calculated)','FontSize',10);
% xl.FontName = defFont{1};
% 
% grid on
% grid minor
% xlim(limH)
% ylim(limH)

%% NEW FIGURE FOR D/H
marksize = 18;
defFont = {'Royal Society of Chemistry'};
tYpos = limH(2)-0.02;

linCol = [166,217,106]./255;
cycCol = [253,174,97]./255;

h = figure;
set(h,'Units','Centimeters','Position',[2 2 8.3 7])

plot(xfitH,yfitHlin,'Color',linCol,'LineWidth',1.5); hold on
plot(xfitH,PI95lin,':','Color',linCol,'LineWidth',1.5); hold on

plot(xfitH,yfitHlin+dyCIlin,'-','Color',linCol,'LineWidth',0.5); hold on
plot(xfitH,yfitHlin-dyCIlin,'-','Color',linCol,'LineWidth',0.5); hold on

set(gca,'FontSize',8)

for j = 1:length(xdatHlin)
    h1 = scatter(xdatHlin(j),ydatHlin(j),marksize,'filled');
    h1.MarkerFaceColor = linCol;
    h1.MarkerEdgeColor = 'k';
end

% if m(2) > 0
%     t = text(limH(1)+0.02,tYpos,sprintf('y = %.3fx + %.3f, R^2 = %.2f',m(1),m(2),rsq));
% elseif m(2) < 0
%     t = text(limH(1)+0.02,tYpos,sprintf('y = %.3fx - %.3f, R^2 = %.2f',m(1),abs(m(2)),rsq));
% end
% t.FontName = defFont{1};
% t.FontSize = 7.5;
% t.Color = 'k';
hold on

f = get(gca,'Children');
f = flipud(f);
% l = legend([f(19) f(4)],{'Linear','Cyclic'});
% l.Location = 'southeast';
% l.Box = 'on';
% l.EdgeColor = 'w';

axLabels = {'0.85','0.90','0.95','1.00','1.05'};
set(gca,'XTickLabel',axLabels,'YTickLabel',axLabels)
xl = xlabel('EIF - Experimental','FontSize',10);
xl.FontName = defFont{1};
yl = ylabel('EIF - Calculated','FontSize',10);
xl.FontName = defFont{1};

text(1.01,0.83,'Linear','FontSize',10)

grid on
grid minor
xlim(limH)
ylim(limH)


ax_2 = axes('Units','centimeters','Position',[1.3 3.9 2.4 2.4]);
axes(ax_2)

plot(xfitH,yfitHcyc,'Color',cycCol,'LineWidth',1.5); hold on
plot(xfitH,PI95cyc,':','Color',cycCol,'LineWidth',1.5); hold on

plot(xfitH,yfitHcyc+dyCIcyc,'-','Color',cycCol,'LineWidth',0.5); hold on
plot(xfitH,yfitHcyc-dyCIcyc,'-','Color',cycCol,'LineWidth',0.5); hold on

for j = 1:length(xdatHcyc)
    h1 = scatter(xdatHcyc(j),ydatHcyc(j),marksize,'filled');
    h1.MarkerFaceColor = cycCol;
    h1.MarkerEdgeColor = 'k';
end

ax_2.XTickLabel = [];
ax_2.YTickLabel = [];

text(0.82,1.02,'Cyclic','FontSize',10)

grid on
% grid minor
xlim(limH)
ylim(limH)

%% D/H plot with no separation of linear and cyclic molecules
limH = [min([min(xdatH) min(ydatH)])-0.03...
        max([max(xdatH) max(ydatH)])+0.03];
xfitH = linspace(limH(1),limH(2),100);

colcode = keto(36:end,1);

[s1H,s2H] = fit(xdatH,ydatH,'poly1');
mH       = [s1H.p1 s1H.p2];
yfitH   = mH(1).*xfitH + mH(2);
rsqH     = s2H.rsquare;
[sa1H,sa2H] = polyfit(xdatH,ydatH,1);
[~,dyCIH] = polyconf(sa1H,xfitH,sa2H,'predopt','curve','Alpha',1-pilevel);

marksize = 18;
defFont = {'Royal Society of Chemistry'};
tYpos = limH(2)-0.02;

HCol(1,:) = [166,217,106]./255;
HCol(2,:) = [253,174,97]./255;

h = figure;
set(h,'Units','Centimeters','Position',[2 2 8.3 7])

plot(xfitH,yfitH,'Color','k','LineWidth',1); hold on
plot(xfitH,yfitH+dyCIH,'--','Color','k','LineWidth',0.5); hold on
plot(xfitH,yfitH-dyCIH,'--','Color','k','LineWidth',0.5); hold on

set(gca,'FontSize',8)

for j = 1:length(xdatH)
    h1 = scatter(xdatH(j),ydatH(j),marksize,'filled');
    if strcmp(colcode(j),'C')
        h1.MarkerFaceColor = HCol(1,:);
    else
        h1.MarkerFaceColor = HCol(2,:);
    end
    h1.MarkerEdgeColor = 'k';
end

t = text(limH(1)+0.02,tYpos,sprintf('y = %.3fx + %.3f, R^2 = %.2f',mH(1),mH(2),rsqH));
t.FontName = defFont{1};
t.FontSize = 7.5;
t.Color = 'k';
hold on

f = get(gca,'Children');
f = flipud(f);
l = legend([f(19) f(4)],{'Linear','Cyclic'});
l.Location = 'southeast';
l.Box = 'on';
l.EdgeColor = 'k';

axLabels = {'0.85','0.90','0.95','1.00','1.05'};
set(gca,'XTickLabel',axLabels,'YTickLabel',axLabels)
xl = xlabel('EIF (Experimental)','FontSize',10);
xl.FontName = defFont{1};
yl = ylabel('EIF (Calculated)','FontSize',10);
xl.FontName = defFont{1};

grid on
xlim(limH)
ylim(limH)
