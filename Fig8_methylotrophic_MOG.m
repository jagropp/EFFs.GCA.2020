function [] = Fig8_methylotrophic_MOG

% Function to calculate isotopic fractionations in methylotrophic
% methanogenesis and to generate Fig. 8 in Gropp et al., GCA, 2020.
%
% This function finds a numerical solution for steady state carbon isotope
% fractionation in the methylotrophic pathway, for the reactions:
% (1) CH3OH/CH3-S-CoM
% (2) CH3-S-CoM/CH4
% (3) CH3-S-CoM/--/CO2
% See details and explanation in appendix B.2 in the paper.

prompt1 = 'Type temperature in degrees Celsius (0-120): ';
prompt2 = 'Type number of simulations (integers only, 1 or larger): ';

% Calculate equilibrium fractionation factors
Tc        = input(prompt1);
beta_vals = calc_132_betas(Tc);
a13eq(1)  = beta_vals.beta13_values(5)/beta_vals.beta13_values(14);
a13eq(2)  = beta_vals.beta13_values(14)/beta_vals.beta13_values(4);
a13eq(3)  = beta_vals.beta13_values(14)/beta_vals.beta13_values(2);

sims      = round(input(prompt2)); % Number of simulations, preferably larger than 500
if sims < 1
    sims = 1;
end
min_a13kf = 0.95; % Minimal C KFF value of reactions 1 and 3
max_a13kf = 0.97; % Maximal C KFF value of reactions 1 and 3
min_f     = -3;   % log10 of the minimal reversibility for reactions 1 and 2
max_f     = 0;    % log10 of the maximal reversibility for reactions 1 and 2
phi_net   = 1e-7; % Net rate (arbitrary and unitless)
R_A       = 1e-6; % Initial isotopic composition of compound A (arbitrary)
Rred_ox   = [0.95 0.75 0.5]; % i.e., 20:1, 3:1 and 1:1

a13kf  = zeros(sims,3);   % Forward KFFs
a13kr  = zeros(sims,3);   % Backward KFFs
f_vals = zeros(sims,3);   % Reversibility along the pathway
R_B    = zeros(sims,3);   % 13R of CH3-S-CoM
R_C    = zeros(sims,3);   % 13R of CH4
R_D    = zeros(sims,3);   % 13R of CO2
kln13a = zeros(2,sims,3); % 1000ln13alpha

reverseStr = '';
for i = 1:sims
    for j = 1:3
        % Draw KFFs and reversibility from defined distributions
        a13kf(i,:)  = (max_a13kf-min_a13kf).*rand(1,3) + min_a13kf;
        a13kf(:,2)  = 0.9615; % Value from Scheller 2013
        a13kr(i,:)  = a13eq.*a13kf(i,:);
        f_vals(i,:) = 10.^((max_f-min_f).*rand(1,3) + min_f);
        f_vals(i,3) = 0.75;
        jf          = phi_net./(1-f_vals(i,:));
        jr          = phi_net.*f_vals(i,:)./(1-f_vals(i,:));
        
        C = [1 1 1].*1e-6';
        options = odeset('NonNegative',1:3,'RelTol',1e-1,'AbsTol',1e-9);
        [~,C] = ode15s(@frac_hydr_ODE,[0 1e10],C,options,jf,jr,a13kf(i,:),a13kr(i,:),R_A,Rred_ox(j));
        R_B(i,j) = C(end,1);
        R_C(i,j) = C(end,2);
        R_D(i,j) = C(end,3);
        kln13a(1,i,j) = 1000*log(R_A./R_C(i,j));
        kln13a(2,i,j) = 1000*log(R_A./R_D(i,j));
    end
    percentDone = 100 * i / sims;
    msg = sprintf('Progress: %3.2f percent\n', percentDone);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end

%% PLOT FIGURE 8 FROM PAPER

xlabel_values = {['1000ln^{13}\alpha_{CH_3OH' char(8211) 'CH_4} (' char(8240) ')'],...
                ['1000ln^{13}\alpha_{CH_3OH' char(8211) 'CO_2} (' char(8240) ')']};
title_values = {'20:1','3:1','1:1'};
figure(1)
set(1,'Units','Centimeters','position',[10 2 22 10])

gapx      = 0.7;
gapy      = 1.9;
gap2      = 0.2;
figWidth  = 10;
figHeight = 2.5;

ha(1) = axes('Units','Centimeters','Position',[gapx gapy+2*gap2+2*figHeight figWidth figHeight]);
ha(2) = axes('Units','Centimeters','Position',[gapx gapy+gap2+figHeight figWidth figHeight]);
ha(3) = axes('Units','Centimeters','Position',[gapx gapy figWidth figHeight]);
ha(4) = axes('Units','Centimeters','Position',[2*gapx+figWidth gapy+2*gap2+2*figHeight figWidth figHeight]);
ha(5) = axes('Units','Centimeters','Position',[2*gapx+figWidth gapy+gap2+figHeight figWidth figHeight]);
ha(6) = axes('Units','Centimeters','Position',[2*gapx+figWidth gapy figWidth figHeight]);

binWidth = 35;
pos = [1 2 3; 4 5 6];
for i = 1:2
    for j = 1:3
        axes(ha(pos(i,j)))
        if i == 1
            binVals = linspace(30,100,binWidth);
        elseif i == 2
            binVals = linspace(-50,30,binWidth);
        end 
        histogram(kln13a(i,:,j),binVals,'LineStyle','none','FaceColor','k')
        hold on
        if j == 3; xlabel(xlabel_values{i}); end
        if j < 3;  set(gca,'XTick',[]); end
        set(gca,'FontSize',13,'YTick',[]);
        if i == 1 && j == 1
            title('R_{r/o} = 20:1',...
            'Units','normalized','Position',[0.84 0.6 0])
        else
            title(title_values{j},...
            'Units','normalized','Position',[0.92 0.7 0])
        end
    end
end

local_path = '/Users/jonathag/Dropbox (Weizmann Institute)/Apps/Overleaf/EFF paper/figures/';
print('-f1',[local_path,'methylotrophic'],'-dpng','-r400');
print('-f1',[local_path,'methylotrophic'],'-depsc');

function ddt = frac_hydr_ODE(~,C,jf,jr,akf,akr,Ra,Rred_ox)
% Funtion for the ODE solver
    
ddt = zeros(size(C));

Rb = C(1);
Rc = C(2);
Rd = C(3);

A = 1;
B = Rred_ox;
C = A-B;

fnet = jf(1)-jr(1);

ba = A*jr(1)*akr(1)*Rb;
bc = B*jf(2)*akf(2)*Rb;
bd = C*jf(3)*akf(3)*Rb;

ab = A*jf(1)*akf(1)*Ra;
cb = B*jr(2)*akr(2)*Rc;
db = C*jr(3)*akr(3)*Rd;

ddt(1) = ab + cb + db - (ba + bc + bd);
ddt(2) = bc - (cb + B*fnet*Rc);
ddt(3) = bd - (db + C*fnet*Rd);
end

end
