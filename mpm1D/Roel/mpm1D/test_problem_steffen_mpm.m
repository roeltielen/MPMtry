% Test problem Steffen

% This file provides input for MPM to compute numerical solution, computes
% an exact solution, illustrates results and provides the RME error.

clear all
close all
beep off
clc 

%% Input
% Constants
density = 1E2;
Youngs_modulus = 1E2;
gravitational_acceleration = 0; 
load = 0;
Length = 1;
tau = 1;
Heigth = 1.15*Length;

% Mesh properties
number_elements = 10; 
number_particles_per_element = 3.5 %3 + 3/(number_elements);  
element_size = Heigth/number_elements; 
mesh = 0:Heigth/(number_elements):Heigth;

% Particles are equally distrubuted over bar.
pos_p_glob = zeros(number_elements*number_particles_per_element,1);
pos_p_loc = zeros(number_elements*number_particles_per_element,1);

for p = 1:number_particles_per_element*number_elements
    pos_p_glob(p) = p/(number_particles_per_element*number_elements);  
end
clear p 

% Adjust the local positions of the particles
for p = 1:number_elements*number_particles_per_element
    pos_p_loc(p) = mod(pos_p_glob(p),element_size)/element_size;
end
clear p

% Time step 
CFL_number = 0.5;
total_time = 4; 
t_cr = element_size/sqrt(Youngs_modulus/density);
t_step = CFL_number*t_cr;
number_time_steps = floor(total_time/t_step); 
t = 0:t_step:(number_time_steps-1)*t_step;

% Initial conditions 
% particles
displacement_p_initial = zeros(number_elements*number_particles_per_element,1);
velocity_p_initial = zeros(number_elements*number_particles_per_element,1);
stress_p_initial = zeros(number_elements*number_particles_per_element,1);
% nodes
displacement_initial = zeros(number_elements+1,1); 
velocity_initial = zeros(number_elements+1,1); 

%% FLAGS
% Boundary conditions
both_ends_fixed = 0;

% Update volume? Yes: 1; No: 0.
volume_update = 1;

% UNDER CONSTRUCTION: Use a lumped matrix or consistent matrix? Lumped: 1; Consistent: 0.
lumped = 1;

% Change global positions? Yes: 1; No: 0.
change_glob_pos = 1;

% Change local positions? Yes: 1; No: 0.
change_loc_pos = 1;

% Identical or modified Lagrangian algorithm? Identical: 1; Modified: 0.
langranian = 1;

% Use momentumvector to determine nodal accelerations? Yes: 1; No: 0.
momentum = 0;

% Update nodal positions? Yes: 1; No: 0.
ULFEM = 0;

% Reset mesh after one timestep? Yes: 1; No: 0.
reset = 1;

% Show a pulse for every grid crossing? Yes: 1; No: 0.
pulse = 1;

% UNDER CONSTRUCTION: Linear or Quadratic shape functions? Linear: 1; Quadratic: 0.
shape = 1;

%% Flags
% Plot displacement versus time for the selected particle? Yes: 1; No: 0 
displ_time = 1; 
% Plot velocity versus time for thr selected particle? Yes: 1; No: 0 
velocity_time = 1;
% Plot velocity vector versus time for the selected particle? Yes: 1; No: 0
velo_time_nodes = 1;


%% Compute the solution using MPM
[displacement_mpm, velocity_mpm, Mass,F_int_plot,stress_particle,indicator,velocity_mpm_nodes] = MPM_1D_Linear_vibrating_string_free_end(density,...
    Youngs_modulus, gravitational_acceleration, load,Heigth, Length,...
    number_elements, element_size, number_particles_per_element,...
    pos_p_glob, pos_p_loc, t_step, number_time_steps, total_time, ...
    displacement_p_initial, velocity_p_initial, stress_p_initial,...
    displacement_initial, velocity_initial, both_ends_fixed, ...
    change_glob_pos,tau,langranian,reset);



%% Plot displacement/velocity versus time for one particle
% Select the particle and node to plot
particle_number = floor(number_elements*number_particles_per_element);
node_number = floor((number_elements+1)/4)

% Obtain exact solution for particles
for p = 1:number_elements*number_particles_per_element
[displacement_exact(p,:) velocity_exact(p,:)]  = exact_solution_dynamic_bar(density,tau,Length,pos_p_glob(p),t,number_time_steps);
end
clear p 

% Obtain exact solution
for n = 1:number_elements + 1
[displacement_exact_nodes(n,:) velocity_exact_nodes(n,:)]  = exact_solution_dynamic_bar(density,tau,Length,mesh(n),t,number_time_steps);
end
clear n

% Plot the displacement of the particle
if displ_time == 1
    figure(1)
    plot_particle_displacement_vs_time(particle_number, displacement_mpm...
        (particle_number,:), displacement_exact(particle_number,:),t, ULFEM, change_loc_pos)
end

% Plot the velocity of the particle
if velocity_time == 1
    figure(2)
    plot_node_velocity_vs_time(particle_number,...
        velocity_mpm(particle_number,:),...
        velocity_exact(particle_number,:),t, ULFEM, change_loc_pos)
end

if velo_time_nodes == 1
       figure(3)
       plot(t,velocity_mpm_nodes(node_number,:),'LineWidth',2)
       hold on
       plot(t,velocity_exact_nodes(node_number,:),'--r','LineWidth',2)
       xlabel('time [s]','Fontsize',12)
       ylabel('velocity [m/s]', 'Fontsize',12)
       title(sprintf('Velocity of node %d',node_number),'FontSize', 12)
       legend('MPM')
end

figure(4)
plot(t,F_int_plot,'b','LineWidth',1.5)
hold on
plot(t,indicator(node_number - 1,:),'r','LineWidth',1.5)
plot(t,indicator(node_number,:),'r','LineWidth',1.5)
plot(t,indicator(node_number + 1,:),'r','LineWidth',1.5)
xlabel('Time [s]','Fontsize',12)
ylabel('Force [N]','Fontsize',12)
legend('Grid crossing','F_{int}')
title(sprintf('Internal force of node 11'),'FontSize', 12)
axis([0 2 -0.6 0.6])


%% Accuracy of MPM solution compared to exact solution at timestep t_check

% Define time to check solution
t_check = floor(length(t)/4);

% Determine RMS error
Error = norm(displacement_exact(:,floor(t_check))-displacement_mpm(:,floor(t_check)))...
    /sqrt(number_elements*number_particles_per_element)

t_step = CFL_number*t_cr



