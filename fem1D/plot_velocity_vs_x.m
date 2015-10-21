function [] = plot_velocity_vs_x(T, num_sol, exact_sol, mesh) 
% Plot the displacement corresponding to a particular node versus time 

% INPUT: time for which the displacement should be plotted, numerical
% solution, exact solution

plot(mesh,num_sol, 'LineWidth',2)
hold on
plot(mesh,exact_sol,'--r', 'LineWidth',2)
xlabel('position [m]', 'FontSize', 12)
set(gca,'FontSize',11)
ylabel('velocity [m/s]','FontSize', 12)
title(sprintf('Velocity at time %g [s]',T)) 
legend('FEM', 'Exact', 'Location','SouthEast')
hold on
set(0,'DefaultFigureColor',[1 1 1])

set(gcf, 'PaperPosition', [0 0 6 5]);
set(gcf, 'PaperSize', [6 5]);
%saveas(gcf,'vel_x_oedom_p0.pdf')