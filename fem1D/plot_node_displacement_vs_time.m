function [] = plot_node_displacement_vs_time(node_number,...
    num_sol, exact_sol, t) 
% Plot the displacement corresponding to a particular node versus time 

% INPUT: nodes for which the displacement should be plotted, numerical
% solution, exact solution, time vector

plot(t,num_sol, 'LineWidth',2)
hold on
plot(t,exact_sol,'--r', 'LineWidth',2)
xlabel('time [s]', 'FontSize', 12)
set(gca,'FontSize',11)
ylabel('displacement [m]','FontSize', 12)
title(sprintf('Displacement of node %d',node_number),'FontSize', 12)
legend('FEM', 'Exact')
hold on
set(0,'DefaultFigureColor',[1 1 1])

set(gcf, 'PaperPosition', [0 0 6 5]);
set(gcf, 'PaperSize', [6 5]);
%saveas(gcf,'displ_time_oedom_p0.pdf')



