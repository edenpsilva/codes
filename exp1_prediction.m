%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright Weifeng Liu  
 
%CNEL
%July 1, 2008
%
%description:
%compare the perfomance of LMS and KLMS for Mackey Glass time series
%one step prediction
%Learning curves
%
%Usage:
%ch2, m-g prediction
%
%Outside functions called:
%none
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all,
close all
clc
%======filter config=======
%time delay (embedding) length
TD = 10;
%kernel parameter
a = 1;%fixed
%noise std
np =.04;
%data size
N_tr = 500;
N_te = 100;%
%%======end of config=======


disp('Learning curves are generating. Please wait...');

%======data formatting===========
load MK30   %MK30 5000*1
ruido = randn(size(MK30));%%ruido de distribuição normal
MK30 = MK30+np*ruido;%adicionando ruido aos dados de entada
MK30 = MK30 - mean(MK30); %centralizando os dados de entrada

%500 training data
train_set = MK30(1501:4500);

%100 testing data
test_set = MK30(4601:4900);

%data embedding
X = zeros(TD,N_tr);
for k=1:N_tr
    X(:,k) = train_set(k:k+TD-1)';
end
T = train_set(TD+1:TD+N_tr);

X_te = zeros(TD,N_te);
for k=1:N_te
    X_te(:,k) = test_set(k:k+TD-1)';
end
T_te = test_set(TD+1:TD+N_te);
%======end of data formatting===========

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
%ploting best aproximation
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
% % colorbar

% colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
linhas = 30;
title('Correlação variando o alfa e a taxa de aprendizagem')
hold on
contourf(lr_ks, alpha, correlation,linhas)
xlabel ('taxa de aprendizagem')
ylabel ('alfa')
colorbar

%saving data
save('exp2_1.mat')