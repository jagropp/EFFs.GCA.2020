% Code to reproduce Fig. S3 in Gropp et al., GCA, 2020 of the relation
% between 12CH2D2 and 13CH3D. We used here the equations presented in the
% paper by Cao et al., 2019, GCA (DOI: 10.1016/J.GCA.2019.01.021), and
% compare our calculations with the predicted EFFs are presented in the
% paper.

% Set number of simulations
p = 1e5;

% Run local functions to calcualte 12CH2D2 and 13CH3D
[x_Cao,y_Cao,k_Cao] = Cao2019(p);
[x,y,k]             = Cao2019_Gropp(p);

% Calculate clumped isotopologue distributions for T = 0-700oC
beta = [ 0.0307,-0.3703,1.8021,-1.4258,0.3210;
        -0.1104, 0.7348,1.1508,-2.8365,1.3194];
Tk = (0:700)+273.15;
temp = beta(:,1).*1e12./Tk.^4 + beta(:,2).*1e9./Tk.^3 + ...
                  beta(:,3).*1e6./Tk.^2 + beta(:,4).*1e3./Tk + beta(:,5);
eq_13CH3D  = temp(1,:);
eq_12CH2D2 = temp(2,:);

% Plot figure
set(figure,'Units','Centimeters','Position',[10 2 15 12])

color_mat = [228,26,28
             55,126,184
             77,175,74]./255;

han_eq = plot(eq_13CH3D,eq_12CH2D2,'k');
hold on

for i = 1:3
    plot(x_Cao(i,k_Cao{i}),y_Cao(i,k_Cao{i}),'Color',color_mat(i,:),...
        'LineStyle',':','LineWidth',2)
    hold on
end
         
for i = 1:3
    han(i) = plot(x(i,k{i}),y(i,k{i}),'Color',color_mat(i,:),'LineWidth',2);
    hold on
end
legend([han_eq han(1) han(2) han(3)],...
    {'Equilibrium','[1,0,0,0]','[1,1,0,0]','[1,1,1,0]'},...
    'Box','off','Location','southeast','FontSize',14)

xlim([-1 6.8])
ylim([-145 24])
xlabel(['\Delta^{13}CH_3D (' char(8240) ')'])
ylabel(['\Delta^{12}CH_2D_2 (' char(8240) ')'])
set(gca,'FontSize',14)

function [x,y,k] = Cao2019_Gropp(p)
% Calculate clumped isotopologue compositions with the calculated EFFs
% presented in Gropp et al., GCA, 2020.

% Calculated EFFs at 25oc
a2eq_CH_F420  = 1.011;
a2eq_CH_CoB   = 2.050;
a2eq_CH2_F420 = 1.1193;
a2eq_CH2_CoB  = 2.2697;
a2eq_CH3_CoB  = 1.9902;

gamma_13D     = random('Uniform',1.000,1.000,1,p);
gamma_DD      = random('Uniform',1.000,1.000,1,p);
gamma_13D_p   = random('Uniform',0.998,1.000,1,p);
gamma_DD_p    = random('Uniform',0.994,1.000,1,p);
KIE_s         = random('Uniform',0.820,1.000,1,p);
KIE_s3        = random('Uniform',0.810,0.900,1,p); % Scheller et al. 2013
KIE_p         = random('Uniform',0.400,1.000,1,p);
KIE_p3        = random('Uniform',0.375,0.450,1,p); % Scheller et al. 2013
Delta_13Deq_1 = 4.5596/1000;
Delta_13Deq_2 = 4.6920/1000;
Delta_13Deq_3 = 5.2191/1000;
Delta_DDeq_2  = 13.3790/1000;
Delta_DDeq_3  = 15.6063/1000;

% [1,0,0,0]
x(1,:) = 1000.*(((gamma_13D.*KIE_s.*(1+Delta_13Deq_1) + ...
    gamma_13D_p.*KIE_p.*(1/a2eq_CH_F420) + gamma_13D_p.*KIE_p.*(1/a2eq_CH_F420) + ...
    gamma_13D_p.*KIE_p3.*(1/a2eq_CH_CoB))./...
    (KIE_s + KIE_p.*(1/a2eq_CH_F420) + KIE_p.*(1/a2eq_CH_F420) + ...
    KIE_p3.*(1/a2eq_CH_CoB))) - 1);

y(1,:) = 1000.*(((gamma_DD_p.*KIE_p.*KIE_s.*(1/a2eq_CH_F420) + ...
    gamma_DD_p.*KIE_p.*(KIE_s + KIE_p.*(1/a2eq_CH_F420)).*(1/a2eq_CH_F420) + ...
    gamma_DD_p.*KIE_p3.*(KIE_s + KIE_p.*(1/a2eq_CH_F420) + ...
    KIE_p.*(1/a2eq_CH_F420)).*(1/a2eq_CH_CoB))./...
    ((6/16).*(KIE_s + KIE_p.*(1/a2eq_CH_F420) + ...
    KIE_p.*(1/a2eq_CH_F420) + KIE_p3.*(1/a2eq_CH_CoB)).^2))-1);

% [1,1,0,0]
x(2,:) = 1000.*(((2*gamma_13D.*KIE_s.*(1+Delta_13Deq_2) + ...
    gamma_13D_p.*KIE_p.*(1/a2eq_CH2_F420) + gamma_13D_p.*KIE_p3.*(1/a2eq_CH2_CoB))./...
    (2.*KIE_s + KIE_p.*(1/a2eq_CH2_F420) + KIE_p3.*(1/a2eq_CH2_CoB)))-1);

y(2,:) = 1000*((gamma_DD.*KIE_s.^2.*(1+Delta_DDeq_2) + ...
    2.*gamma_DD_p.*KIE_s.*KIE_p.*(1/a2eq_CH2_F420) ...
    + gamma_DD_p.*KIE_p3.*(2.*KIE_s + KIE_p.*(1/a2eq_CH2_F420)).*(1/a2eq_CH2_CoB))./...
    ((6/16).*((2.*KIE_s + KIE_p.*(1/a2eq_CH2_F420) + KIE_p3.*(1/a2eq_CH2_CoB)).^2))-1);

% [1,1,1,0]
x(3,:) = 1000.*((3*gamma_13D.*KIE_s3.*(1+Delta_13Deq_3)+...
    gamma_13D_p.*KIE_p3./a2eq_CH3_CoB)./...
    (3.*KIE_s3 + KIE_p3./a2eq_CH3_CoB) - 1);

y(3,:) = 1000.*((8.*KIE_s3.*(gamma_DD.*KIE_s3.*...
    (1+Delta_DDeq_3)+gamma_DD_p.*KIE_p3./a2eq_CH3_CoB))./...
    (((3.*KIE_s3 + KIE_p3./a2eq_CH3_CoB).^2)) - 1);

for i = 1:3
    k{i} = convhull([x(i,:)',y(i,:)']);
end

end 

function [x,y,k] = Cao2019(p)
% Produce results as presented in Cao et al., GCA, 2019. Equations and
% parameters as in Tables 1 and 2 there.

gamma_13D     = random('Uniform',1.000,1.000,1,p);
gamma_DD      = random('Uniform',1.000,1.000,1,p);
gamma_13D_p   = random('Uniform',0.998,1.000,1,p);
gamma_DD_p    = random('Uniform',0.994,1.000,1,p);
D_KIE         = random('Uniform',0.820,1.000,1,p);
Delta_13Deq_i = random('Uniform',0.005,0.005,1,p);
Delta_DDeq_i  = random('Uniform',0.018,0.018,1,p);
R_1           = random('Uniform',0.300,0.800,1,p);

% [1,0,0,0]
x(1,:) = 1000.*(((gamma_13D.*D_KIE.*(1+Delta_13Deq_i) + 3.*gamma_13D_p.*R_1)./...
    (D_KIE + 3.*R_1)) - 1);
y(1,:) = 1000.*(((8.*gamma_DD_p.*R_1.*(D_KIE + R_1))./((D_KIE + 3.*R_1).^2)) - 1);

% [1,1,0,0]
x(2,:) = 1000.*((gamma_13D.*D_KIE.*(1+Delta_13Deq_i)+gamma_13D_p.*R_1)./...
    (D_KIE + R_1) - 1);

y(2,:) = 1000.*(((2.*(D_KIE.*(gamma_DD.*D_KIE.*(1+Delta_DDeq_i)+...
    4.*gamma_DD_p.*R_1)+gamma_DD_p.*(R_1.^2)))./...
    (3.*(D_KIE + R_1).^2)) - 1);
    
% [1,1,1,0]
x(3,:) = 1000.*((3*gamma_13D.*D_KIE.*(1+Delta_13Deq_i)+gamma_13D_p.*R_1)./...
    (3.*D_KIE + R_1) - 1);

y(3,:) = 1000.*(1.*((8.*D_KIE.*(gamma_DD.*D_KIE.*(1+Delta_DDeq_i)+gamma_DD_p.*R_1))./...
    (((3.*D_KIE + R_1).^2)) - 1));

for i = 1:3
    k{i} = convhull([x(i,:)',y(i,:)']);
end

end


