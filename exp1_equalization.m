%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright Aaron Liu 
%CNEL
%July 1, 2008
%
%description:
%compare the performance of LMS and KLMS in nonlinear channel equalization
%learning curve
%
%Usage:
%ch2, nonlinear channel equalization, figure 2-5
%
%Outside functions called:
%None
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all,
close all
clc
tic
%======filter config=======
%time delay (embedding) length
TD = 5;
D = 2;
h = .1;%kernel size
%======end of config=======

%=========data===============
% Generate binary data
u = randn(1,2500)>0;
u = 2*u-1;

% Nonlinear channel
z = u+0.5*[0,u(1:end-1)];
% Channel noise
ns = 0.4*randn(1,length(u));
% Ouput of the nonlinear channel
y = z - 0.9*z.^2 + ns;

%data size
N_tr = 1000;
N_te = 50;

%data embedding
X = zeros(TD,N_tr);
for k=1:N_tr
    X(:,k) = y(k:k+TD-1)';
end
% Test data
X_te = zeros(TD,N_te);
for k=1:N_te
    X_te(:,k) = y(k+N_tr:k+TD-1+N_tr)';
end

% Desired signal
T = zeros(N_tr,1);
for ii=1:N_tr
    T(ii) = u(D+ii);
end

T_te = zeros(N_te,1);
for ii=1:N_te
    T_te(ii) = u(D+ii+N_tr);
end
%======end of data===========

%=======init================
mse_te = zeros(1,N_tr);
mse_te_k = zeros(1,N_tr);
mse_te_ks = zeros(1,N_tr);
%=======end of init=========


%=========KSIG===================

%learning rate (adjustable)
alpha = 0.025:0.025:2;
lr_ks = 0.01:0.01:1;

na = length(alpha);%size(alfa,1);%length(alpha)
nl = length(lr_ks);%size(lr_ks,1);

%init
e_ks = zeros(N_tr,1);
y = zeros(N_tr,1);
mse_te_ks = zeros(N_tr,1);

% n=1 init
e_ks(1) = T(1);
y(1) = 0;
mse_te_ks(1) = mean(T_te.^2);
aux = zeros (N_tr-1,1);%auxliar varable, it's used to calculatin less
correlation = zeros(na,nl);%correlation vector
%%monitoring execution time
tic
total = na*nl;
atual = 0;
clc
fprintf('porcentagem->%3.2f%%\n',atual/total*100);
fprintf('tempo: %1.1f minutos\n',toc/60);

for i=1:na
    for j = 1:nl
        atual = atual + 1;    
        %%training
        eta = lr_ks(j)*alpha(i);%new eta
        aux(1) = tanh(alpha(i)*e_ks(1));%new coefficienT
        for n=2:N_tr

            %training
            ii = 1:n-1;
            y(n) = eta*aux(ii)'*(exp(-sum((X(:,n)*ones(1,n-1)-X(:,ii)).^2)))';
            e_ks(n) = T(n) - y(n);
            aux(n) = tanh(alpha(i)*e_ks(n));%new coefficient
        end
        y_te = zeros(N_te,1);
        for jj = 1:N_te
            y_te(jj) = eta*(aux(1:n))'*(exp(-sum((X_te(:,jj)*ones(1,n)-X(:,1:n)).^2)))';
        end
        corr = corrcoef(T_te, y_te);
        correlation(i,j) = corr(1,2);
    end 
        clc
        fprintf('porcentagem->%1.2f%%\n',atual/total*100);
        fprintf('tempo: %1.1f minutos\n',toc/60);

end



maximo = max(correlation(:));
[inda, indl] = find(correlation==maximo);%find max correlation

fprintf('máxima correlacao: %1.4f\n',maximo);
fprintf('alpha: %1.1f; learning rate:  %1.2f\n',alpha(inda), lr_ks(indl));
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%with the best values of alpha and lr, plot to see the graphic differences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%training to find the value to the best aproximation
eta = lr_ks(indl)*alpha(inda);
for n=2:N_tr
    aux(n-1) = tanh(alpha(inda)*e_ks(n-1));
    %training
    ii = 1:n-1;
    y(n) = eta*aux(ii)'*(exp(-sum((X(:,n)*ones(1,n-1)-X(:,ii)).^2)))';
    e_ks(n) = T(n) - y(n);
end
y_te = zeros(N_te,1);
for jj = 1:N_te
    y_te(jj) = eta*(e_ks(1:n))'*(exp(-sum((X_te(:,jj)*ones(1,n)-X(:,1:n)).^2)))';
end


% %ploting best aproximation
% figure
% plot(y_te,':k')
% hold on
% plot(T_te,'b')
% legend('aproximação', 'sinal teste')
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% surf(lr_ks, alpha, correlation)
% xlabel ('taxa de aprendizagem')
% ylabel ('alfa')
% zlabel ('correlação')
% colorbar

figure
linhas = 30;
set(gca, 'FontSize', 14);
set(gca, 'FontName', 'Arial');
title('correlation by alpha versus learning rate')
hold on
contourf(lr_ks, alpha, correlation,linhas)
colorbar
colormap gray

xlabel ('learning rate')
ylabel ('alpha')


%saving data
save('exp1_equalization')

