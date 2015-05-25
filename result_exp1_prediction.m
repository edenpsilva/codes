%cleanning 
clc
clear

load 'exp1_prediction'

%   mse_te_k = zeros(1,N_tr);
%   h=1;
% %=========Kernel LMS===================
%     lr_k = .2;
%     %init
%     e_k = zeros(N_tr,1);
%     y = zeros(N_tr,1);
%     y_te = zeros(N_te,1);
%     % n=1 init
%     e_k(1) = T(1);
%     y(1) = 0;
%     mse_te_k(1) = mean(T_te.^2);
%     % start
%     for n=2:N_tr
%         %training
%         ii = 1:n-1;
%         y(n) = lr_k*e_k(ii)'*(exp(-sum((X(:,n)*ones(1,n-1)-X(:,ii)).^2)*h))';
% 
%         e_k(n) = T(n) - y(n);
% 
%         %testing
%         y_te = zeros(N_te,1);
%         for jj = 1:N_te
%             ii = 1:n;
%             y_te(jj) = lr_k*e_k(ii)'*(exp(-sum((X_te(:,jj)*ones(1,n)-X(:,ii)).^2)*h))';
%         end
%         err = T_te - y_te;
%         mse_te_k(n) = mean(err.^2);
%        
%     end
%     corrk = corrcoef(T_te, y_te);
% 
%     %=========end of Kernel LMS================
% 
% 
% 
% 
% maximo = max(correlation(:));
% [inda, indl] = find(correlation==maximo);%find max correlation
% 
% fprintf('máxima correlacao: %1.4f\n',maximo);
% fprintf('alpha: %1.1f; learning rate:  %1.2f\n',alpha(inda), lr_ks(indl));
%     
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %with the best values of alpha and lr, plot to see the graphic differences
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%training to find the value to the best aproximation
% eta = lr_ks(indl)*alpha(inda);
% for n=2:N_tr
%     aux(n-1) = tanh(alpha(inda)*e_ks(n-1));
%     %training
%     ii = 1:n-1;
%     y(n) = eta*aux(ii)'*(exp(-sum((X(:,n)*ones(1,n-1)-X(:,ii)).^2)))';
%     e_ks(n) = T(n) - y(n);
% end
% y_te = zeros(N_te,1);
% for jj = 1:N_te
%     y_te(jj) = eta*(e_ks(1:n))'*(exp(-sum((X_te(:,jj)*ones(1,n)-X(:,1:n)).^2)))';
% end
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
% xlabel ('learning rate')
% ylabel ('alpha')
% zlabel ('correlation')
% % colorbar

% colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
