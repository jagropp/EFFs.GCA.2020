% 11.04.2020
%
% Produce Fig. 2 from Cao et al., 2019.

clear
p = 1e5;

gamma_13D     = random('Uniform',1.000,1.000,1,p);
gamma_DD      = random('Uniform',1.000,1.000,1,p);
gamma_13D_p   = random('Uniform',0.998,1.000,1,p);
gamma_DD_p    = random('Uniform',0.994,1.000,1,p);
D_KIE         = random('Uniform',0.820,1.000,1,p);
Delta_13Deq_i = random('Uniform',0.005,0.005,1,p);
Delta_DDeq_i  = random('Uniform',0.018,0.018,1,p);
R_1           = random('Uniform',0.300,0.800,1,p);

% 1000
x(1,:) = 1000.*(((gamma_13D.*D_KIE.*(1+Delta_13Deq_i) + 3.*gamma_13D_p.*R_1)./...
    (D_KIE + 3.*R_1)) - 1);
y(1,:) = 1000.*(((8.*gamma_DD_p.*R_1.*(D_KIE + R_1))./((D_KIE + 3.*R_1).^2)) - 1);

% 1100
x(2,:) = 1000.*((gamma_13D.*D_KIE.*(1+Delta_13Deq_i)+gamma_13D_p.*R_1)./...
    (D_KIE + R_1) - 1);

y(2,:) = 1000.*(((2.*(D_KIE.*(gamma_DD.*D_KIE.*(1+Delta_DDeq_i)+...
    4.*gamma_DD_p.*R_1)+gamma_DD_p.*(R_1.^2)))./...
    (3.*(D_KIE + R_1).^2)) - 1);
    
% 1110
x(3,:) = 1000.*((3*gamma_13D.*D_KIE.*(1+Delta_13Deq_i)+gamma_13D_p.*R_1)./...
    (3.*D_KIE + R_1) - 1);

y(3,:) = 1000.*(1.*((8.*D_KIE.*(gamma_DD.*D_KIE.*(1+Delta_DDeq_i)+gamma_DD_p.*R_1))./...
    (((3.*D_KIE + R_1).^2)) - 1));

% nbins = 30;
% Xedges = linspace(-1,5,nbins+1);
% Yedges = linspace(-140,10,nbins+1);
% for i = 1:3
%     [N(:,:,i)] = histcounts2(x(i,:),y(i,:),Xedges,Yedges);
%     Z(:,:,i) = N(:,:,i)'./p;
% end
% Z(Z==0) = NaN;
% x_vals = movmean(Xedges(2:end),2);
% y_vals = movmean(Yedges(2:end),2);
% [X,Y] = meshgrid(x_vals,y_vals);

for i = 1:3
    k{i} = convhull([x(i,:)',y(i,:)']);
end

% figure
% subplot(1,2,1)
for i = 1:3
%     scatter(x(i,:),y(i,:),'filled'); hold on
    plot(x(i,k{i}),y(i,k{i}),'Color',color_mat(i,:))
    hold on
end
xlim([-1 5])
ylim([-145 10])
xlabel('\Delta^{13}CH_3D')
ylabel('\Delta^{12}CH_2D_2')

