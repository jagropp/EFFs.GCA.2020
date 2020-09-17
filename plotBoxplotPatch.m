clear
load('190514_exp_calc_CNOH_dat')
% dat1(:,5) = NaN;
dat1(:,:) = dat1([6 5 3 1 2 4 7 10 8 9 11],:);
dat2(:,:) = dat2([6 5 3 1 2 4 7 10 8 9 11],:);

c = [178,24,43
214,96,77
244,165,130
253,219,199
247,247,247
209,229,240
146,197,222
67,147,195
33,102,172]./255;

ci(1:7,:) = repmat(c(7,:),7,1);
ci(8:11,:) = repmat(c(7,:),4,1);

cii(1:7,:) = repmat(c(2,:),7,1);
cii(8:11,:) = repmat(c(2,:),4,1);

wbox = 0.3;
xA(:,1) = linspace(0.6,11.6,12);
xA(:,2) = xA(:,1)+wbox;
xB(:,1) = linspace(1.1,12.1,12);
xB(:,2) = xB(:,1)+wbox;

prc1i = prctile(dat1,[25 50 75],2);
prc1IQR = prc1i(:,3)-prc1i(:,1);
prc1highol = prc1i(:,3)+prc1IQR.*1.5;
prc1lowol  = prc1i(:,1)-prc1IQR.*1.5;
cond1i = find(dat1 >= prc1highol | dat1 <= prc1lowol);
dat1i = dat1;
prc1 = prctile(dat1i,[25 50 75],2);

prc2i = prctile(dat2,[25 50 75],2);
prc2IQR = prc2i(:,3)-prc2i(:,1);
prc2highol = prc2i(:,3)+prc2IQR.*1.5;
prc2lowol  = prc2i(:,1)-prc2IQR.*1.5;
cond2i = find(dat2 >= prc1highol | dat2 <= prc2lowol);
dat2i = dat2;
prc2 = prctile(dat2i,[25 50 75],2);

h = figure;
set(h,'Units','Centimeters','Position',[2 2 17.1 9])

for i = 1:11
    yyaxis left
    cond1a = find(dat1(i,:) <= prc1highol(i) & dat1(i,:) >= prc1lowol(i));
    cond1b = find(dat1(i,:) >= prc1highol(i) | dat1(i,:) <= prc1lowol(i));
    prc1(i,:) = prctile(dat1(i,cond1a),[25 50 75],2);
    line([mean([xA(i,1) xA(i,2)]) mean([xA(i,1) xA(i,2)])],...
        [min(dat1(i,cond1a)) max(dat1(i,cond1a))],'Color','k','LineWidth',0.5)
    hold on
    patch([xA(i,:) flip(xA(i,:))],[prc1(i,1) prc1(i,1) prc1(i,3) prc1(i,3)],ci(i,:))
    hold on
    line(xA(i,:),[prc1(i,2) prc1(i,2)],'Color','k','LineWidth',0.5)
    hold on
    line(xA(i,:),[nanmean(dat1(i,:)) nanmean(dat1(i,:))],'Color','k','LineStyle',':')
    hold on
    scatter(ones(1,length(dat1(i,cond1b)))*mean([xA(i,1) xA(i,2)]),dat1(i,cond1b),...
        'Marker','+')
    ylabel(['CNO absolute deviation [' char(8240) ']'],'FontName','Acumin pro')
    ylim([0 18])
    
    yyaxis right
    cond2a = find(dat2(i,:) <= prc2highol(i) & dat2(i,:) >= prc2lowol(i));
    cond2b = find(dat2(i,:) >= prc2highol(i) | dat2(i,:) <= prc2lowol(i));
    prc2(i,:) = prctile(dat2(i,cond2a),[25 50 75],2);
    line([mean([xB(i,1) xB(i,2)]) mean([xB(i,1) xB(i,2)])],...
        [min(dat2(i,cond2a)) max(dat2(i,cond2a))],'Color','k','LineWidth',0.5)
    hold on
    patch([xB(i,:) flip(xB(i,:))],[prc2(i,1) prc2(i,1) prc2(i,3) prc2(i,3)],cii(i,:))
    hold on
    line(xB(i,:),[prc2(i,2) prc2(i,2)],'Color','k','LineWidth',0.5)
    hold on
    line(xB(i,:),[nanmean(dat2(i,:)) nanmean(dat2(i,:))],'Color','k','LineStyle',':')
    hold on
    scatter(ones(1,length(dat2(i,cond2b)))*mean([xB(i,1) xB(i,2)]),dat2(i,cond2b),...
        'Marker','+')
    ylabel(['H absolute deviation [' char(8240) ']'],'FontName','Acumin pro')
end
% line([7.5 7.5],[0 180],'Color',[1 1 1].*0.5,'LineWidth',0.5)
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10 11])
set(gca,'XTickLabel',{'\tauHCTH_{D3BJ}','\tauHCTH','HCTH','M06-L_{D3}',...
    'M06-L','OLYP','SOGGA11','M11-L','O3LYP','\tauHCTH_{hyb}','B3LYP'})
set(gca,'FontSize',9)
xtickangle(45)
xlim([0 12])



