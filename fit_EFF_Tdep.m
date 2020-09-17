% b = [1;1;1;1];

x = linspace(0,700,351)+273.15;
% x = linspace(0,374,188)+273.15;

modelfun = @(b,x)(b(1).*1e12./x.^4 + b(2).*1e9./x.^3 + ...
                  b(3).*1e6./x.^2 + b(4).*1e3./x + b(5));
beta0 = [1;1;1;1;1];
clear b beta
for i = 1:size(y,1)
    beta(i,:) = nlinfit(x,y(i,:),modelfun,beta0);
end

y1 = beta(:,1).*1e12./x.^4 + beta(:,2).*1e9./x.^3 + ...
                  beta(:,3).*1e6./x.^2 + beta(:,4).*1e3./x + beta(:,5);
              
% figure;
% plot(1e6./x.^2,y,'o',1e6./x.^2,y1,'--')
% set(gca,'FontSize',14)

% b(:,1) = beta(:,1);
% b(:,2) = beta(:,2);
% b(:,3) = beta(:,3);
% b(:,4) = beta(:,4);

% for C
b(:,1) = beta(:,1).*1e6;
b(:,2) = beta(:,2).*1e5;
b(:,3) = beta(:,3).*1e4;
b(:,4) = beta(:,4).*1e4;

% % For H
% b(:,1) = beta(:,1).*1e3;
% b(:,2) = beta(:,2).*1e2;
% b(:,3) = beta(:,3).*1e2;
% b(:,4) = beta(:,4).*1e2;

% % For CH
% b(:,1) = beta(:,1).*1e3;
% b(:,2) = beta(:,2).*1e2;
% b(:,3) = beta(:,3).*1e2;
% b(:,4) = beta(:,4).*1e2;

% % For HH
% b(:,1) = beta(:,1).*1e1;
% b(:,2) = beta(:,2).*1e0;
% b(:,3) = beta(:,3).*1e0;
% b(:,4) = beta(:,4).*1e0;

b(:,5) = beta(:,5);





