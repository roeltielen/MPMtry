function [] = plot_particle_displacement_vs_time(particle_number,...
    num_sol, exact_sol, t, ULFEM, change_loc_pos) 
% Plot the displacement corresponding to a particular node versus time 

% INPUT: nodes for which the displacement should be plotted, numerical
% solution, exact solution, time vector

plot(t,num_sol, 'LineWidth',2)
hold on
plot(t,exact_sol,'--r', 'LineWidth',2)
xlabel('time [s]', 'FontSize', 12)
set(gca,'FontSize',11)
ylabel('displacement [m]','FontSize', 12)
title(sprintf('Position of particle %d',particle_number),'FontSize', 12)
if ULFEM == 1
        legend('ULFEM','Exact')
    else if change_loc_pos == 0
            legend('FEM','Exact')
         else
            legend('MPM','Exact')
         end
    end
hold on
set(0,'DefaultFigureColor',[1 1 1])

set(gcf, 'PaperPosition', [0 0 6 5]);
set(gcf, 'PaperSize', [6 5]);
saveas(gcf,'displ_time_mpm.pdf')