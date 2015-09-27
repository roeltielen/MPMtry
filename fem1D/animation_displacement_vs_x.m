function [] = animation_displacement_vs_x(T, num_sol, exact_sol, mesh,...
    max_displ, min_displ, height) 
% Plot the displacement corresponding to a particular node versus time 

% INPUT: time for which the displacement should be plotted, numerical
% solution, exact solution

plot(mesh,num_sol, 'LineWidth',2)
hold on
plot(mesh,exact_sol,'--r', 'LineWidth',2)
ylim([min_displ, max_displ])
xlim([-0.1, 1.1*height])
xlabel('position [m]', 'FontSize', 12)
set(gca,'FontSize',11)
ylabel('displacement [m]','FontSize', 12)
title(sprintf('Displacement at time %f [s]',T)) 
legend('FEM', 'Exact')
hold on
set(0,'DefaultFigureColor',[1 1 1])

set(gcf, 'PaperPosition', [0 0 6 6]);
set(gcf, 'PaperSize', [6 6]);