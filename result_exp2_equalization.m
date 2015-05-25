clear all
clc
load exp2_equalization

%fprintf('Ploting MSE\n');
% ==========plot and test===================
%%%%%%%%%
% figure
% % plot(mse_te,'b-','LineWidth',2)
% % hold on
% plot(mse_te_k,'r--','LineWidth',2)
% hold on
% plot(mse_te_ks,'k:','LineWidth',2)
% set(gca, 'FontSize', 14);
% set(gca, 'FontName', 'Arial');
% legend('KLMS', 'KSIG')
% xlabel('iteration')
% ylabel('testing MSE')

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