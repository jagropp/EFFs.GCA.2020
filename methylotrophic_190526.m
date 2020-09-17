addpath('/Users/jonathag/Dropbox (Weizmann Institute)/Jonathan_(with_Itay)/Code/general_functions/')
clear
% methylottrophic pathway - last update 11.02.2020

% (1)   CH3OH/CH3-S-CoM
% (2)   CH3-S-CoM/CH4
% (3-8) CH3-S-CoM/--/CO2

beta_vals = calc_132_betas(25);
aeq(1) = beta_vals.beta13_values(5)/beta_vals.beta13_values(14);
aeq(2) = beta_vals.beta13_values(14)/beta_vals.beta13_values(4);
aeq(3) = beta_vals.beta13_values(14)/beta_vals.beta13_values(13);
aeq(4) = beta_vals.beta13_values(13)/beta_vals.beta13_values(11);
aeq(5) = beta_vals.beta13_values(11)/beta_vals.beta13_values(10);
aeq(6) = beta_vals.beta13_values(10)/beta_vals.beta13_values(9);
aeq(7) = beta_vals.beta13_values(9)/beta_vals.beta13_values(8);
aeq(8) = beta_vals.beta13_values(8)/beta_vals.beta13_values(2);
aeq(9) = 1.0;

p      = 5e1;
minakf = 0.95;
maxakf = 0.97;
minf   = 1e-4;
maxf   = 1;
phinet = 1e-2;
Ra     = 1e-8;
Rred_ox = [0.95 0.75 0.6];

akf  = zeros(p,9);
akr  = zeros(p,9);
f    = zeros(p,9);
Rb   = zeros(p,3); % CH3-S-CoM
Rc   = zeros(p,3);
Rd   = zeros(p,3);
Re   = zeros(p,3);
Rf   = zeros(p,3);
Rg   = zeros(p,3);
Rh   = zeros(p,3);
Ri   = zeros(p,3);
Rj   = zeros(p,3);
klna = zeros(3,p,3);

for i = 1:p
    for j = 1:3
        akf(i,:) = (maxakf-minakf).*rand(1,9) + minakf;
        akf(:,2) = 1/1.04;
        akr(i,:) = aeq.*akf(i,:);
        f(i,:) = (maxf-minf).*rand(1,9) + minf;
%         f(i,1) = 1e-5; % Reversibility of Mta
%         f(i,2) = 1e-5; % Reversibility of Mcr
        f(i,3:9) = ones(1,7).*0.999;
%         f(i,9) = 1e-6;
        jf = phinet./(1-f(i,:));
        jr = phinet.*f(i,:)./(1-f(i,:));
        
        C = [1 1 1 1 1 1 1 1 1].*1e-6';
        options = odeset('NonNegative',1:9,'RelTol',1e-1,'AbsTol',1e-9);
        [t,C] = ode15s(@frac_hydr_ODE,[0 1e10],C,options,jf,jr,akf(i,:),akr(i,:),Ra,Rred_ox(j));
        Rb(i,j) = C(end,1);
        Rc(i,j) = C(end,2);
        Rd(i,j) = C(end,3);
        Re(i,j) = C(end,4);
        Rf(i,j) = C(end,5);
        Rg(i,j) = C(end,6);
        Rh(i,j) = C(end,7);
        Ri(i,j) = C(end,8);
        Rj(i,j) = C(end,9);
        klna(i,1,j) = 1000*log(Ra./Rc(i,j));
        klna(i,2,j) = 1000*log(Ra./Ri(i,j));
        klna(i,3,j) = 1000*log(Ra./Rj(i,j));
        fprintf('%d out of %d\n',i,p)
    end
end

title_values = {['ln\alpha(CH_3OH/CH_4) (' char(8240) ')'],...
                ['ln\alpha(CH_3OH/CO_2) (' char(8240) ')'],...
                ['ln\alpha(CH_3OH/Biomass) (' char(8240) ')']};

q   = 100;
gap = 3;
FaceAlphaVal = 0.5;
colors_code = [228,26,28; 55,126,184; 77,175,74]./255;
set(figure,'Units','Centimeters','position',[10 2 28 9])
ha = tight_subplot(1, 3, [0.05 0.05], [0.18 0.02], 0.02);

for i = 1:3
    for j = 1:3
        pd(i,j) = fitdist(klna(:,i,j),'Kernel');
        x(i,:,j) = linspace(min(pd(i,j).InputData.data)-gap,...
            max(pd(i,j).InputData.data)+gap,q);
        y(i,:,j) = pdf(pd(i,j),x(i,:,j));
%         subplot(1,3,i)
        axes(ha(i))
%         hold on
%         histogram(klna(:,i,j),20,'LineStyle','none','Normalization','probability',...
%                 'FaceAlpha',FaceAlphaVal,'FaceColor',colors_code(j,:))
        hold on
        plot(x(i,:,j),y(i,:,j),'Color',colors_code(j,:),'LineWidth',2)
        xlabel(title_values{i})
        set(gca,'FontSize',13,'YTick',[])
    end
end

l = legend(num2str(Rred_ox(1)),num2str(Rred_ox(2)),num2str(Rred_ox(3)));
l.Location = 'northeast';
l.Box      = 'off';
title(l,'R_{red/ox}')

function ddt = frac_hydr_ODE(~,C,jf,jr,akf,akr,Ra,Rred_ox)
ddt = zeros(size(C));

Rb = C(1);
Rc = C(2);
Rd = C(3);
Re = C(4);
Rf = C(5);
Rg = C(6);
Rh = C(7);
Ri = C(8);
Rj = C(9);

% Rred_ox = 0.5;               % Default 0.75;
Biomass_sink = 0.05;
A = 1;                        % Default 4
B = (1-Biomass_sink)*Rred_ox; % Default 3
C = A-B-Biomass_sink;         % Default 1

fnet = jf(1)-jr(1);

ab = A*jf(1)*akf(1)*Ra;
ba = A*jr(1)*akr(1)*Rb;
bc = B*jf(2)*akf(2)*Rb;
cb = B*jr(2)*akr(2)*Rc;
bd = C*jf(3)*akf(3)*Rb;
db = C*jr(3)*akr(3)*Rd;
de = C*jf(4)*akf(4)*Rb;
ed = C*jr(4)*akr(4)*Re;
ef = C*jf(5)*akf(5)*Re;
fe = C*jr(5)*akr(5)*Rf;
fg = C*jf(6)*akf(6)*Rf;
gf = C*jr(6)*akr(6)*Rg;
gh = C*jf(7)*akf(7)*Rg;
hg = C*jr(7)*akr(7)*Rh;
hi = C*jf(8)*akf(8)*Rh;
ih = C*jr(8)*akr(8)*Ri;

bj = Biomass_sink*jf(9)*akf(9)*Rb;
jb = Biomass_sink*jr(9)*akr(9)*Rj;

ddt(1) = ab + cb + db + jb - (ba + bc + bd + bj);
ddt(2) = bc - (cb + B*fnet*Rc);
ddt(3) = bd + ed - (db + de);
ddt(4) = de + fe - (ed + ef);
ddt(5) = ef + gf - (fe + fg);
ddt(6) = fg + hg - (gf + gh);
ddt(7) = gh + ih - (hg + hi);
ddt(8) = hi - (ih + C*fnet*Ri);
ddt(9) = bj - (jb + Biomass_sink*fnet*Rj);

end
