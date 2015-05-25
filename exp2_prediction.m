%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Experimento 2 - teste 2 -  com alpha fixo, 
% tax de aprendizagem dado por 
%mu_sig = mu_lms*2*sperança(sech(alfa*ruido).^2);xd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright Weifeng Liu  
 
%CNEL
%July 1, 2008
%
%description:
%compare the performance of LMS and KLMS for Mackey Glass time series
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

tic
disp('Learning curves are generating. Please wait...');
mc = 1;%test number
%correlation4 = zeros(4,mc);%store the correlation

alpha = 0.025:0.025:2;
learning_rates_KSIG= zeros(1,mc);
correlation = zeros(2, mc);
erro_klms = zeros(mc, N_tr);
erro_ksig = zeros(mc, N_tr);
for i=1:mc %in every iteratiom, the data change the noise
    %======data formatting===========
    load MK30   %MK30 5000*1
    ruido = np*randn(size(MK30));%%ruido de distribuição normal
    MK30 = MK30+ruido;%adicionando ruido aos dados de entada
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

     
    %=========Kernel LMS===================

    %learning rate (adjustable)
    %    lr_k = .1;
    lr_k = .2;
    %   lr_k = .6;

    %init
    e_k = zeros(N_tr,1);
    y = zeros(N_tr,1);
    mse_te_k = zeros(N_tr,1);

    % n=1 init
    e_k(1) = T(1);
    y(1) = 0;
    mse_te_k(1) = mean(T_te.^2);
    % start
    for n=2:N_tr
        %training
        ii = 1:n-1;
        y(n) = lr_k*e_k(ii)'*(exp(-sum((X(:,n)*ones(1,n-1)-X(:,ii)).^2)))';
        e_k(n) = T(n) - y(n);

        %testing MSE
        y_te = zeros(N_te,1);
        for jj = 1:N_te
            y_te(jj) = lr_k*e_k(1:n)'*(exp(-sum((X_te(:,jj)*ones(1,n)-X(:,1:n)).^2)))';
        end
        err = T_te - y_te;
        %erro_klms(i,n) = err;
        mse_te_k(n) = mean(err.^2);
        erro_klms(i,n) = mse_te_k(n);

    end
    corrk = corrcoef(T_te, y_te);
    correlation(1,i)=corrk(1,2);
%     
    %=========end of Kernel LMS================

    %=========KSIG===================

    %learning rate 
    inda = 35;
    faixaruido = ruido((1500 + n):(1500 + TD + n -1));
    lr_ks = lr_k*2*(mean((sech(faixaruido*alpha(inda))).^2)*mean(faixaruido.^2))/mean((tanh(faixaruido*alpha(inda))).^2);
    learning_rates_KSIG(i) = lr_ks;
    %finding the best alpha
    clc
    toc
    fprintf('%d de %d\n',i,mc)
    y_te_ks = zeros(N_te,1);
    %init
    e_ks = zeros(N_tr,1);
    y_ks = zeros(N_tr,1);
    mse_te_ks = zeros(N_tr,1);

    % n=1 init
    e_ks(1) = T(1);
    y_ks(1) = 0;
    mse_te_ks(1) = mean(T_te.^2);
    aux = zeros (N_tr-1,1);
    aprox =10;
    % start
    eta = lr_ks*alpha(inda);
    aux(1) = tanh(alpha(inda)*e_ks(1));
    for n=2:N_tr
        
        %training
        ii = 1:n-1;
        %eta é lr_ks*alpha(inda)
        y_ks(n) = eta*aux(ii)'*(exp(-sum((X(:,n)*ones(1,n-1)-X(:,ii)).^2)))';

        e_ks(n) = T(n) - y_ks(n);
        aux(n) = tanh(alpha(inda)*e_ks(n));

        %testing MSE
        y_te_ks = zeros(N_te,1);
        for jj = 1:N_te
            y_te_ks(jj) = eta*(aux(1:n))'*(exp(-sum((X_te(:,jj)*ones(1,n)-X(:,1:n)).^2)))';
        end
         err = T_te - y_te_ks;
         %erro_klms(i,n) = err;
        mse_te_ks(n) = mean(err.^2);
        erro_ksig(i,n) = mse_te_ks(n);
    end
    corrks = corrcoef(T_te, y_te_ks);
    correlation(2,i)= corrks(1,2);
 
    %=========end of Kernel sigmoide================

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ploting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% plot(mse_te_k,'r--','LineWidth',2);
% 
% hold on
% plot(mse_te_ks,'b:','LineWidth',2);
% 
% 
% set(gca, 'FontSize', 14);
% set(gca, 'FontName', 'Arial');
% 
% legend('KLMS', 'KSIG')
% xlabel('iteração')
% ylabel('MSE')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure
% hold on
% plot(correlation(1,:),'r--','LineWidth',2);
% 
% hold on
% plot(correlation(2,:),'b:','LineWidth',2);
% 
% 
% set(gca, 'FontSize', 14);
% set(gca, 'FontName', 'Arial');
% 
% legend('KLMS', 'KSIG')
% xlabel('iteração')
% ylabel('correlação')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure
% plot(learning_rates_KSIG,'LineWidth', 2)
% xlabel('iteração')
% ylabel('sem unidade definida')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure
% 
% % semilogy(std(erro_klms),'k:')
% % hold on
% % semilogy(std(erro_ksig), 'r')
% % grid
% 
%  semilogy(mean(erro_klms),'k:','LineWidth', 2)
%  hold on
%  semilogy(mean(erro_ksig), 'r','LineWidth', 2)
%  grid
% xlabel('iteração')
% ylabel('média dos mse')
% legend('KLMS', 'KSIG')
figure

plot(mean(erro_klms(:,2:end)),'k','LineWidth', 2)
hold on
plot(mean(erro_ksig(:,2:end)), 'k-.','LineWidth', 2)
grid

set(gca, 'FontSize', 14);
set(gca, 'FontName', 'Arial');
xlabel('iteration')
ylabel('MSE')
legend('KLMS', 'KSIG')
colormap gray

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% fprintf('Testes: %d\n',mc)
% fprintf('---------------------------------------\n')
% fprintf('Média da melhor correlação:\n')
% fprintf('KLMS - %1.4f\n',mean(correlation(1,:)))
% fprintf('KSIG - %1.4f\n',mean(correlation(2,:)))
% fprintf('---------------------------------------\n')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% fprintf('alfa - %1.4f\n',alpha(inda))
% fprintf('---------------------------------------\n')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% maximo = max(correlation(2,:));
% indmax = find(correlation(2,:)==maximo);
% minimo = min(correlation(2,:));
% indmin= find(correlation(2,:)==minimo);
% 
% fprintf('máxima correlação: %1.4f, lr: %1.4f\n',maximo,learning_rates_KSIG(indmax))
% fprintf('mínima correlação: %1.4f, lr: %1.4f\n',minimo,learning_rates_KSIG(indmin))
% fprintf('---------------------------------------\n')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% maximo = max(learning_rates_KSIG);
% indmax = find(learning_rates_KSIG==maximo);
% minimo = min(learning_rates_KSIG);
% indmin= find(learning_rates_KSIG==minimo);
% 
% fprintf('correlação: %1.4f, máxima taxa de aprendizagem: %1.4f\n',correlation(2,indmax),learning_rates_KSIG(indmax))
% fprintf('correlação: %1.4f, mínimo taxa de aprendizagem: %1.4f\n',correlation(2,indmin),learning_rates_KSIG(indmin))
% fprintf('---------------------------------------\n')
% 
% 
 save('exp2_prediction.mat')