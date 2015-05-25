%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                   
%Experimento 2 - teste 2 -  com alpha fixo, 
% tax de aprendizagem dado por 
%mu_sig = mu_lms*2*speran�a(sech(alfa*ruido).^2);xd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all

load exp2_prediction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
% xlabel('itera��o')
% ylabel('MSE')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure
% hold on
% plot(correlation(1,:),'b--','LineWidth',2);
% 
% hold on
% plot(correlation(2,:),'k-','LineWidth',2);
% 
% 
% set(gca, 'FontSize', 14);
% set(gca, 'FontName', 'Arial');
% 
% legend('KLMS', 'KSIG')
% xlabel('itera��o')
% ylabel('correla��o')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure
% plot(learning_rates_KSIG,'LineWidth', 2)
% xlabel('itera��o')
% ylabel('sem unidade definida')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure

% semilogy(std(erro_klms),'k:','LineWidth', 2)
% hold on
% semilogy(std(erro_ksig), 'r','LineWidth', 2)
% grid

%  semilogy(mean(erro_klms),'k:','LineWidth', 2)
%  hold on
%  semilogy(mean(erro_ksig), 'r','LineWidth', 2)
%  grid
% xlabel('itera��o')
% ylabel('m�dia dos mse')
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% fprintf('Testes: %d\n',mc)
% fprintf('---------------------------------------\n')
% fprintf('M�dia da melhor correla��o:\n')
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
% fprintf('m�xima correla��o: %1.4f, lr: %1.4f\n',maximo,learning_rates_KSIG(indmax))
% fprintf('m�nima correla��o: %1.4f, lr: %1.4f\n',minimo,learning_rates_KSIG(indmin))
% fprintf('---------------------------------------\n')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% maximo = max(learning_rates_KSIG);
% indmax = find(learning_rates_KSIG==maximo);
% minimo = min(learning_rates_KSIG);
% indmin= find(learning_rates_KSIG==minimo);
% 
% fprintf('correla��o: %1.4f, m�xima taxa de aprendizagem: %1.4f\n',correlation(2,indmax),learning_rates_KSIG(indmax))
% fprintf('correla��o: %1.4f, m�nimo taxa de aprendizagem: %1.4f\n',correlation(2,indmin),learning_rates_KSIG(indmin))
% fprintf('---------------------------------------\n')
