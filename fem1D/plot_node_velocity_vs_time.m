function [] = plot_node_velocity_vs_time(node_number,num_sol, exact_sol, t) 
% Plot the displacement corresponding to a particular node versus time 

% INPUT: nodes for which the displacement should be plotted, numerical
% solution, exact solution, time vector

plot(t,num_sol, 'LineWidth',2)
hold on
plot(t,exact_sol,'--r', 'LineWidth',2)
xlabel('time [s]', 'FontSize', 12)
set(gca,'FontSize',11)
ylabel('velocity [m/s]','FontSize', 12)
title(sprintf('Velocity of node %d',node_number)) 
legend('FEM', 'Exact')
hold on
set(0,'DefaultFigureColor',[1 1 1])

set(gcf, 'PaperPosition', [0 0 6 6]);
set(gcf, 'PaperSize', [6 6]);
%saveas(gcf,'vel_time_vs.pdf')