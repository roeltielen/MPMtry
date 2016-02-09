function [ ] = wave_propagation_2D() 
tic
% Input file containing
% 1) material parameters: porosity, density of the water phase, density of 
% the solid phase, Youngs_modulus, bulk_modulus, Darcy permeability, total
% stress, initial pore pressure, height/length of the sample;
% 2) mesh properties: number of elements
% 3) time discretization properties: total time, CFL number
% 4) initial conditions: initial velocity of water, initial velocity of
% solid skeleton, initial total stress, initial pore pressure, initial
% effective stress
% 5) boundary conditions

clear all
close all
beep off
clc  

%% Input
% Material parameters
porosity = 0.4;
density_w = 1E3;
density_s = 2.6E3;
Youngs_modulus = 2.5E9; 
bulk_modulus = 2E9; 
permeability = 1E-5;
total_stress = -10E3; 
pore_pressure = -10E3; 
height = 2;
width = 1;
eff_stress = total_stress - pore_pressure;
grav_accel = 10;

% Mesh properies
n_ey = 2; 
n_ex = 1;

% Time discretization properties
CFL_number = 0.9;
total_time = 3E-3;

% Initial conditions 
velocity_s_initial = zeros(n_ey + 1,1);
velocity_w_initial = zeros(n_ey + 1,1);
total_stress_initial = zeros(n_ey,1);
total_stress_initial(end) = total_stress;
pore_pressure_initial = zeros(n_ey,1);
%pore_pressure_initial(end) = pore_pressure;
eff_stress_initial = zeros(n_ey,1);
%eff_stress_initial(end) = eff_stress;

% Boundary conditions
v_w_bottom = 0;
v_s_bottom = 0;
water_traction = pore_pressure;

%% Perform FEM computation
% Call FEM_2phase
%[pore_pressure, effective_stress] = 
FEM_2phase(porosity,...
    density_w, density_s, Youngs_modulus, bulk_modulus, permeability,...
    total_stress, pore_pressure, eff_stress, height, width, grav_accel,...
    n_ey, n_ex, CFL_number, total_time, velocity_s_initial,...
    velocity_w_initial, total_stress_initial, pore_pressure_initial,...
    eff_stress_initial, v_w_bottom, v_s_bottom, water_traction);
toc

% %% Illustrate results
% % analytical solution for stress
% [t1, AS, AP] = as_wave_propagation(permeability);
% 
% density_sat = (1-porosity)*density_s + porosity*density_w;
% undrained_constr = Youngs_modulus + bulk_modulus/porosity;
% vel_undrained = sqrt(undrained_constr/density_sat);
% damped_factor = sqrt((porosity*Youngs_modulus/bulk_modulus)/...
%     (1-porosity+porosity*Youngs_modulus/bulk_modulus));
% vel_damped = damped_factor*sqrt(bulk_modulus/density_w);
% element_size = height/n_e;
% t_cr = min(element_size/vel_undrained,element_size/vel_damped);
% t_step = 1E-6; %CFL_number*t_cr;
% number_time_steps = floor(total_time/t_step); 
% t = 0:t_step:(number_time_steps-1)*t_step;
% 
% % Generate element connections 
% xvec  = 0:element_size:height; %global mesh
% elements = zeros(n_e,2);
% for i = 1:n_e
%    elements(i,:) = [xvec(i) xvec(i+1)]; 
% end
% 

% figure(2);
% plot(t,effective_stress(floor(3*n_ey/4),:)/total_stress,'LineWidth',2)
% elements(floor(3*n_ey/4),:)
% hold on
% plot(t1,AS,'--k','LineWidth',2)
% xlabel('time [s]', 'FontSize', 12)
% set(gca,'FontSize',11)
% ylabel('normalized effective stress [-]','FontSize', 12)
% %title(sprintf('Displacement at time %g [s]',T)) 
% legend('FEM', 'Exact', 'Location','southeast')
% % hold on
% set(0,'DefaultFigureColor',[1 1 1])
% set(gcf, 'PaperPosition', [0 0 6 6]);
% set(gcf, 'PaperSize', [6 6]);
% 
% 
% figure(3);
% plot(t,pore_pressure(floor(3*n_e/4),:)/total_stress,'r','LineWidth',2)
% hold on
% plot(t1,AP,'--k','LineWidth',2)
% xlabel('time [s]', 'FontSize', 12)
% set(gca,'FontSize',11)
% ylabel('normalized pore presssure [-]','FontSize', 12)
% % title(sprintf('Displacement at time %g [s]',T)) 
% legend('FEM', 'Exact', 'Location','southeast')
% % hold on
% set(0,'DefaultFigureColor',[1 1 1])
% set(gcf, 'PaperPosition', [0 0 6 6]);
% set(gcf, 'PaperSize', [6 6]);
% 
% % figure(4);
% % x = element_size/2:element_size:n_e*element_size;
% % for nt = 1:20:number_time_steps-1
% %    plot(x',pp(:,nt)) 
% %    pause(1)
% % end


