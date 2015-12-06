clear all
clc
close all

[velocity_mpm, node_number, t] = oedometer_1phase(0,1,1,1,1,0,1,0,1,1,...
    0,1,10);

velocity_ulfem = oedometer_1phase(0,0,1,1,0,1,0,1,0,0,0,1,1);

%grey = [0.4,0.4,0.4];

figure(8)
plot(t,velocity_mpm,'k','LineWidth',2)
hold on
plot(t,velocity_ulfem,'--m','LineWidth',2)
xlabel('time [s]','Fontsize',12)
ylabel('velocity [m/s]', 'Fontsize',12)
title(sprintf('Position of node %d',node_number),'FontSize', 12, 'FontWeight', 'normal')
legend('MPM','ULFEM')
hold on
set(0,'DefaultFigureColor',[1 1 1])
set(gcf, 'PaperPosition', [0 0 16 15]);
set(gcf, 'PaperSize', [16 13]);
