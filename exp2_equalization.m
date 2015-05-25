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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===========testes======
testes = 1;

corrlm = zeros(1, testes);
corrkm = zeros(1, testes);
correlation = zeros(1, testes);
lr_kst = zeros(1, testes);
alpha = 0.025:0.025:2;
inda = 24;
% erro_klms = zeros(1, N_tr);
% erro_ksig = zeros(1, N_tr);

for i = 1:testes
    clc
    fprintf('%d de %d\n',i,testes)
    toc
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
erro_klms = zeros(N_te, N_tr);
erro_ksig = zeros(N_te, N_tr);

 
    %=========Linear LMS===================
    %learning rate (step size)
    lr = .01;%learning rate
    w1 = zeros(1,TD);
    e_l = zeros(N_tr,1);
    for n=1:N_tr
        y = w1*X(:,n);
        e_l(n) = T(n) - y;
        w1 = w1 + lr*e_l(n)*X(:,n)';
        %testing
        if (i == testes)
            err = T_te'-(w1*X_te);
            mse_te(n) = mean(err.^2);
        end
    end
    corrl = corrcoef(T_te,w1*X_te);
    if (isnan(corrl(1,2)))
                corrlm(i) = 0;
    else
                corrlm(i) = corrl(1,2);
    end        
    
    %=========end of Linear LMS================

    %=========Kernel LMS===================
    lr_k = .2;
    %init
    e_k = zeros(N_tr,1);
    y = zeros(N_tr,1);
    y_te = zeros(N_te,1);
    % n=1 init
    e_k(1) = T(1);
    y(1) = 0;
    mse_te_k(1) = mean(T_te.^2);
    % start
    for n=2:N_tr
        %training
        ii = 1:n-1;
        y(n) = lr_k*e_k(ii)'*(exp(-sum((X(:,n)*ones(1,n-1)-X(:,ii)).^2)*h))';

        e_k(n) = T(n) - y(n);

        %testing
        if (i == testes)
            y_te = zeros(N_te,1);
            for jj = 1:N_te
                ii = 1:n;
                y_te(jj) = lr_k*e_k(ii)'*(exp(-sum((X_te(:,jj)*ones(1,n)-X(:,ii)).^2)*h))';
            end

            err = T_te - y_te;
            mse_te_k(n) = mean(err.^2);
            erro_klms(:,n) = err;
        end
    end
    if (i < testes)
        y_te = zeros(N_te,1);
        for jj = 1:N_te
           ii = 1:n;
           y_te(jj) = lr_k*e_k(ii)'*(exp(-sum((X_te(:,jj)*ones(1,n)-X(:,ii)).^2)*h))';
        end
    end
    corrk = corrcoef(T_te,y_te);
     if (isnan(corrk(1,2)))
                corrkm(i) = 0;
    else
                corrkm(i) = corrk(1,2);
    end 


    %=========end of Kernel LMS================

    %find the best learning rate for a fix alpha
    %=========Kernel Sigmoide===================

    %init
    e_ks = zeros(N_tr,1);
    y_ks = zeros(N_tr,1);
    % n=1 init
    e_ks(1) = T(1);
    y_ks(1) = 0;
    mse_te_ks(1) = mean(T_te.^2);
    faixaruido = ns;

   % start
        
    lr_ks = lr_k*2*(mean((sech(faixaruido*alpha(inda))).^2)*mean(faixaruido.^2))/mean((tanh(faixaruido*alpha(inda))).^2);
    %store learning rates
    lr_kst(i) = lr_ks;
    eta = alpha(inda)*lr_ks;
    aux = zeros (N_tr-1,1);
    aprox =10;
    aux(1) = tanh(alpha(inda)*e_ks(1)); 
    for n=2:N_tr
        %training
        ii = 1:n-1;
        y_ks(n) = eta*aux(ii)'*(exp(-sum((X(:,n)*ones(1,n-1)-X(:,ii)).^2)*h))';
        e_ks(n) = T(n) - y_ks(n);
        aux(n) = tanh(alpha(inda)*e_ks(n));
        
        if (i==testes)
            y_te_ks = zeros(N_te,1);
            for jj = 1:N_te
                y_te_ks(jj) = eta*(aux(1:n))'*(exp(-sum((X_te(:,jj)*ones(1,n)-X(:,1:n)).^2)*h))';
            end
             %testing
            err = T_te - y_te_ks;
            mse_te_ks(n) = mean(err.^2);    
            erro_ksig(:,n) = err;
        end
        
    end
    
    if (i<testes)
        y_te_ks = zeros(N_te,1);
        for jj = 1:N_te
            y_te_ks(jj) = eta*(aux(1:n))'*(exp(-sum((X_te(:,jj)*ones(1,n)-X(:,1:n)).^2)*h))';
        end
         %testing
        
    end
    corrks = corrcoef(T_te, y_te_ks);
    if (isnan(corrks(1,2)))
        correlation(i) = 0;
    else
        correlation(i) = corrks(1,2);
    end        
   
end
%end
%=========end of Kernel Sigmoide================

% corr = correlation(:);
% maximo = max(correlation(:));%
% [ind] = find(corr==maximo);%find max correlation
% fprintf('---------------------------------------------------\n')
% fprintf('máxima correlacao: %1.4f\n',maximo);
% fprintf('alpha: %1.2f; learning rate:  %1.2f\n',alpha(inda),lr_kst(ind));
% fprintf('---------------------------------------------------\n')
% fprintf('Correlação média:\n');
% fprintf('LMS  - %1.4f\n',sum(corrlm)/testes);
% fprintf('KLMS - %1.4f\n',sum(corrkm)/testes);
% fprintf('KSIG - %1.4f\n',sum(correlation)/testes);
% fprintf('---------------------------------------------------\n')
% %%%%calculating aproximation to the beste learning rate and alpha
% 
% 
% % ==========plot and test===================
% % figure
% % plot(corrlm,'b-','LineWidth',2);
% % hold on 
% % 
% % plot(corrkm,'r--','LineWidth',2);
% % hold on 
% % plot(correlation,'k:','LineWidth',2);
% % legend('KLMS', 'KSIG')
% % xlabel('iteration')
% % ylabel('correlação')
% % %%%%%%%%%
% figure
% 
% plot(mse_te_k,'r--','LineWidth',2)
% hold on
% plot(mse_te_ks,'k:','LineWidth',2)
% set(gca, 'FontSize', 14);
% set(gca, 'FontName', 'Arial');
% 
% 
% legend('LMS','KLMS', 'KSIG')
% xlabel('iteration')
% ylabel('testing MSE')
% 
% %%%%%%%%%%%%%%%%%
% % % figure
% % % plot(correlation, 'LineWidth',2)
% % % xlabel('teste')
% % % ylabel('correlação')
% 
% %%%%%%%%%%%%%%%%
% figure
% plot(lr_kst,'r','LineWidth',2)
% xlabel('iterações')
% ylabel('taxa de aprendizagem')
% 

%%%%%%%%%
figure
semilogy(std(erro_klms),'k','LineWidth',2)
hold on
semilogy(std(erro_ksig),'k-.','LineWidth',2)
set(gca, 'FontSize', 14);
set(gca, 'FontName', 'Arial');

colormap gray
hold on 
grid

legend('KLMS', 'KSIG')
xlabel('iteration')
ylabel('error')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%saving data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 save('exp2_equalization')

