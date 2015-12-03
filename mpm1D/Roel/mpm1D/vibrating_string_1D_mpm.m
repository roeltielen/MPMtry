% Vibrating string with fixed ends

% This file provides input for MPM to compute numerical solution, computes
% an exact solution, illustrates results and provides the RME error.

clear all
close all
beep off
clc 
format long

%% Input
% Constants
density = 1;
Youngs_modulus = 100;
gravitational_acceleration = 0; 
load = 0;
height = 25; 
alpha = 0;
v_0 = 0.1;

%% FLAGS
% Boundary conditions
both_ends_fixed = 0;

% Update volume? Yes: 1; No: 0.
volume_update = 0;

% UNDER CONSTRUCTION: Use a lumped matrix or consistent matrix? Lumped: 1; Consistent: 0.
lumped = 1;

% Change global positions? Yes: 1; No: 0.
change_glob_pos = 0;

% Change local positions? Yes: 1; No: 0.
change_loc_pos = 0;

% Identical or modified Lagrangian algorithm? Identical: 1; Modified: 0.
lagranian = 1;

% Use momentumvector to determine nodal accelerations? Yes: 1; No: 0.
momentum = 0;

% Update nodal positions? Yes: 1; No: 0.
ULFEM = 0;

% Show a pulse for every grid crossing? Yes: 1; No: 0.
pulse = 0;

% UNDER CONSTRUCTION: Linear or Quadratic shape functions? Linear: 1; Quadratic: 0.
shape = 1;

reset = 0;

deformation = 0;

%% Mesh and particle properties
% Mesh properties
number_elements = 128; 
number_particles_per_element = 4;  
element_size = height/number_elements; 
mesh_one = 0:element_size:height;
mesh_two = 0:element_size/2:height;

% Particles are equally distributed over element.
% initial local position of elements
pos_p_loc = ones(number_elements*number_particles_per_element,1)/(number_particles_per_element+1); 
for i = 1:number_particles_per_element
   pos_p_loc(i:number_particles_per_element:number_elements*number_particles_per_element,1)...
       = i/(number_particles_per_element+1);
end

% initial global position of elements 
pos_p_minloc = zeros(number_elements*number_particles_per_element,1);
for i = 1:number_elements-1 
    pos_p_minloc(i*number_particles_per_element+1:(i+1)*number_particles_per_element,1)...
        = pos_p_minloc(i*number_particles_per_element+1:(i+1)...
        *number_particles_per_element,1) + i*element_size;
end
pos_p_glob = pos_p_minloc + element_size*pos_p_loc; 

% Time step 
CFL_number = 0.9;
total_time = 2.5; 
t_cr = element_size/sqrt(Youngs_modulus/density);
t_step = 1E-3; %CFL_number*t_cr
number_time_steps = floor(total_time/t_step); 
t = 0:t_step:(number_time_steps-1)*t_step;

% Initial conditions 
% particles
displacement_p_initial = zeros(number_elements*number_particles_per_element,1);
if both_ends_fixed == 1
    velocity_p_initial = v_0*sin(pi*pos_p_glob/height);
else
    velocity_p_initial = v_0*sin(pi*pos_p_glob/(2*height));
end
stress_p_initial = zeros(number_elements*number_particles_per_element,1);

% nodes
if shape == 1
    displacement_n_initial = zeros(number_elements + 1,number_time_steps-1);
    if both_ends_fixed == 1
        velocity_n_initial = v_0*sin(pi*mesh_one'/height);
    else
        velocity_n_initial = v_0*sin(pi*mesh_one'/(2*height));
    end
else
    displacement_n_initial = zeros(2*number_elements + 1,number_time_steps-1);
    if both_ends_fixed == 1
        velocity_n_initial = v_0*sin(pi*mesh_one'/height);
    else
        velocity_n_initial = v_0*sin(pi*mesh_one'/(2*height));
    end
end

% Energy vectors
E_kin = zeros(1,number_time_steps);
E_pot = zeros(1,number_time_steps);
E_grav = zeros(1,number_time_steps);
E_trac = zeros(1,number_time_steps);



%% Compute the solution using MPM
if shape == 1
    [displacement_mpm, velocity_mpm,velocity_mpm_nodes, M_lump,...
        displacement_mpm_particles, E_kin,E_pot,E_grav,E_trac,...
        F_diff,F_int_plot,F_ext_plot, stress_particle,indicator] =...
        MPM_1D_Linear_Final(density,...
        Youngs_modulus, gravitational_acceleration, load, height,...
        number_elements, element_size, number_particles_per_element,...
        pos_p_glob, pos_p_loc, t_step, number_time_steps, total_time, ...
        displacement_p_initial, velocity_p_initial, stress_p_initial,...
        displacement_n_initial, velocity_n_initial, both_ends_fixed, ...
        change_glob_pos,change_loc_pos,E_kin,E_pot,E_grav,E_trac,lagranian...
        ,volume_update,lumped, ULFEM, momentum, alpha, reset, deformation);
else
    [displacement_mpm, velocity_mpm,velocity_mpm_nodes, M_lump,...
        displacement_mpm_particles, E_kin,E_pot,E_grav,E_trac,...
        F_diff,F_int_plot,F_ext_plot, stress_particle,indicator] = ...
        MPM_1D_Quadratic_Try(density,...
        Youngs_modulus, gravitational_acceleration, load, height,...
        number_elements, element_size, number_particles_per_element,...
        pos_p_glob, pos_p_loc, t_step, number_time_steps, total_time, ...
        displacement_p_initial, velocity_p_initial, stress_p_initial,...
        displacement_n_initial, velocity_n_initial, both_ends_fixed, ...
        change_glob_pos,change_loc_pos,E_kin,E_pot,E_grav,E_trac,lagranian...
        ,volume_update,lumped, ULFEM, momentum, alpha, reset, deformation);
end

% Initialise position of nodes and particles
position_mpm = zeros(size(displacement_mpm));
position_mpm_particles = zeros(size(displacement_mpm_particles));

% Determine position of nodes and particles
for n = 1:number_time_steps
    position_mpm(:,n) = displacement_mpm(:,n) + mesh_one';
    position_mpm_particles(:,n) = displacement_mpm_particles(:,n) + pos_p_glob;
end
clear n

if both_ends_fixed == 1
    velocity_mpm(:,1) =  v_0*sin(pi*pos_p_glob/height);
    velocity_mpm_nodes(:,1) =  v_0*sin(pi*mesh_one'/height);
else
    velocity_mpm(:,1) =  v_0*sin(pi*pos_p_glob/(2*height));
    velocity_mpm_nodes(:,1) = v_0*sin(pi*mesh_one'/(2*height)) ;
end
% Determine constants
if both_ends_fixed == 1
    w1 = pi*sqrt(Youngs_modulus/density)/height;
    b1 = pi/height;
else
    w1 = pi*sqrt(Youngs_modulus/density)/(2*height);
    b1 = pi/(2*height);
end

%% Obtain the exact solution for the particles(!)
position_exact = zeros(number_elements*number_particles_per_element, number_time_steps);
disp_exact = zeros(number_elements*number_particles_per_element, number_time_steps);
vel_exact = zeros(number_elements*number_particles_per_element, number_time_steps);

if both_ends_fixed == 1
   vel_exact(:,1) = v_0*sin(pi*pos_p_glob/height);
else
   vel_exact(:,1) = v_0*sin(pi*pos_p_glob/(2*height)); 
end

for n = 1:number_time_steps-1
    for p = 1:number_elements*number_particles_per_element
        vel_exact(p,n+1) = v_0*cos(w1*t_step*n)*sin(b1*pos_p_glob(p));
        disp_exact(p,n+1) = v_0/w1*sin(w1*t_step*n)*sin(b1*pos_p_glob(p));    
    end
end

clear node

for n = 1:number_time_steps
    position_exact(:,n) = disp_exact(:,n) + pos_p_glob;
end
clear n

%% Obtain the exact solution for the nodes
position_exact_nodes = zeros(number_elements+1, number_time_steps);
disp_exact_nodes = zeros(number_elements+1, number_time_steps);
vel_exact_nodes = zeros(number_elements+1, number_time_steps);

if both_ends_fixed == 1
   vel_exact_nodes(:,1) = v_0*sin(pi*mesh_one'/height);
else
   vel_exact_nodes(:,1) = v_0*sin(pi*mesh_one'/(2*height)); 
end

for n = 1:number_time_steps-1
    for node = 1:number_elements + 1
        vel_exact_nodes(node,n+1) = 0.1*cos(w1*t_step*n)*sin(b1*mesh_one(node));
        disp_exact_nodes(node,n+1) = 0.1/w1*sin(w1*t_step*n)*sin(b1*mesh_one(node));
    end
end
clear node


for n = 1:number_time_steps
    position_exact_nodes(:,n) = disp_exact_nodes(:,n) + mesh_one';
end
clear n

%% Flags
% Plot displacement versus time for the selected particle? Yes: 1; No: 0 
displ_time = 1; 

% Plot velocity vector versus time for the selected particle? Yes: 1; No: 0
velo_time = 1;

% Plot velocity vector versus time for the selected node? Yes: 1; No: 0
velo_time_nodes = 1;

% Plot energy of the system? Yes: 1; No: 0
energy_plot = 0;

% Plot the internal force of a single node? Yes: 1; No: 0
force_internal_plot = 1;

% Make a waterfall plot of the velocity of all particles? Yes: 1; No: 0
waterfall_plot = 0;

% Make a contour plot of the position of all particles? Yes: 1; No:0
contour_plot = 0;



%% Plot displacement/velocity versus time for one particle
% Select the particle and node
particle_number = floor(number_elements*number_particles_per_element/2);
node_number = floor((number_elements+1)/2);


if displ_time == 1
    figure(1)
    plot_particle_displacement_vs_time(particle_number, position_mpm_particles(particle_number,:), position_exact(particle_number,:),t, ULFEM, change_loc_pos)
    
    figure(8)
    plot(t,position_mpm(node_number+1,:),'b','LineWidth',2)
    hold on
    plot(t,position_exact_nodes(node_number+1,:),'--r','LineWidth',2)
    xlabel('time [s]','Fontsize',12)
    ylabel('position [m]', 'Fontsize',12)
    title(sprintf('Position of node %d',node_number+1),'FontSize', 12)
    if ULFEM == 1
        legend('ULFEM','Exact')
    else if change_loc_pos == 0
            legend('FEM','Exact')
         else
            legend('MPM','Exact')
         end
    end
end

if velo_time == 1
    figure(2)
    plot(t,velocity_mpm(particle_number,:),'b','LineWidth',2)
    hold on
    plot(t,vel_exact(particle_number,:),'--r','LineWidth',2)
    xlabel('time [s]','Fontsize',12)
    ylabel('velocity [m/s]', 'Fontsize',12)
    title(sprintf('Velocity of particle %d',particle_number),'FontSize', 12)
    if ULFEM == 1
        legend('ULFEM','Exact')
    else if change_loc_pos == 0
            legend('FEM','Exact')
         else
            legend('MPM','Exact')
         end
    end
end

if velo_time_nodes == 1
    figure(3)
    plot(t,velocity_mpm_nodes(node_number+1,:),'b','LineWidth',2)
    hold on
    plot(t,vel_exact_nodes(node_number+1,:),'--r','LineWidth',2)
    xlabel('time [s]','Fontsize',12)
    ylabel('velocity [m/s]', 'Fontsize',12)
    title(sprintf('Velocity of node %d',node_number+1),'FontSize', 12)
    if ULFEM == 1
        legend('ULFEM','Exact')
    else if change_loc_pos == 0
            legend('FEM','Exact')
         else
            legend('MPM','Exact')
         end
    end
end

if energy_plot == 1 
    E_total =  E_kin - E_pot + E_trac + E_grav;
    figure(4)
    plot(t,E_kin,'-b','LineWidth',2) 
    hold on
    plot(t,-E_pot,'-r','LineWidth',2) 
    plot(t,E_grav,'-g','LineWidth',2) 
    plot(t,E_trac,'-k','LineWidth',2) 
    plot(t,E_total,'-m','LineWidth',2)
    xlabel('Time [s]','Fontsize',12)
    ylabel('Energy','Fontsize',12)
    legend('E_{kin}','E_{pot}' ,'E_{grav}', 'E_{trac}', 'E_{total}')
end

if force_internal_plot == 1
    figure(5)
    plot(t,F_int_plot,'b','LineWidth',1.5)
    if pulse == 1
        hold on
        plot(t,indicator(node_number,:)/1E5,'m','LineWidth',1.5)
        hold on
        plot(t,indicator(node_number-1,:)/1E5,'r','LineWidth',1.5)
        plot(t,indicator(node_number+1,:)/1E5,'r','LineWidth',1.5)
    end
    xlabel('Time [s]','Fontsize',12)
    ylabel('Force [N]','Fontsize',12)
    legend('F_{int}','Grid crossing','Location','Northwest')
    title(sprintf('Internal force of node %d',node_number),'FontSize', 12)
end

if waterfall_plot == 1
    figure(6)
    waterfall(velocity_mpm)
    hold on
    waterfall(vel_exact)
    xlabel('Time[s]','Fontsize',12)
    ylabel('Particle number','Fontsize',12)
    zlabel('Velocity [m/s]','Fontsize',12)
    title('Velocity of particles over time','Fontsize',12)
    if ULFEM == 1
        legend('ULFEM','Exact')
    else if change_loc_pos == 0
            legend('FEM','Exact')
         else
            legend('MPM','Exact')
         end
    end
end

if contour_plot == 1
    number_contour_plots = 10;
    figure(7)  
    contour(t,pos_p_glob,position_mpm_particles,number_contour_plots)
    hold on
    contour(t,pos_p_glob,position_exact,number_contour_plots)
    xlabel('Time[s]','Fontsize',12)
    ylabel('Position [m]','Fontsize',12)
    title('Position of particles over time','Fontsize',12)
    if ULFEM == 1
        legend('ULFEM','Exact')
    else if change_loc_pos == 0
            legend('FEM','Exact')
         else
            legend('MPM','Exact')
         end
    end
end

%% Accuracy of MPM solution compared to exact solution at timestep t_check

% Define time to check solution
t_check = floor(length(t)/5);

% Determine RME error
Error = norm(position_exact(:,floor(t_check))-position_mpm_particles...
    (:,floor(t_check)))/sqrt(number_elements*number_particles_per_element)







