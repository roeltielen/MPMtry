% Vibrating string with fixed ends

% This file provides input for MPM to compute numerical solution, computes
% an exact solution, illustrates results and provides the global error

clear all
close all
beep off
clc 

%% Input
% Constants
density = 1;
Youngs_modulus = 100;
gravitational_acceleration = 0; 
load = 0;
height = 25; %height/length

% Mesh properties
number_elements = 16; % number of elements
number_particles_per_element = 1; % number of particles per element 
element_size = height/number_elements; % size of an element
pos_p_loc = ones(number_elements*number_particles_per_element,1)*0.5; % initial local positions of particles
mesh = 0:element_size:height;
pos_p_glob = (mesh(1:end-1))' + element_size*pos_p_loc; % initial global position
%pos_p_glob = 1/2*(element_size/number_particles_per_element):...
    %element_size/number_particles_per_element:...
    %height-1/2*(element_size/number_particles_per_element); % initial global position:
% NEEDED FOR INITIAL VELOCITY. If changed, local positions within MPM have 
% to be changed too! 


% Time step 
CFL_number = 0.1;
total_time = 50; 
t_cr = element_size/sqrt(Youngs_modulus/density);
t_step = CFL_number*t_cr;
number_time_steps = floor(total_time/t_step); % set here the total time
t = 0:t_step:(number_time_steps-1)*t_step;

% Initial conditions 
% particles
displacement_p_initial = zeros(number_elements*number_particles_per_element,1);
velocity_p_initial = 0.1*sin(pi/height*pos_p_glob);
stress_p_initial = zeros(number_elements*number_particles_per_element,1);
% nodes
displacement_n_initial = zeros(number_elements+1,1); %initial nodal velocity
velocity_n_initial = zeros(number_elements+1,1); %initial nodal displacement

% Boundary conditions
both_ends_fixed = 1;

% Change global positions? Yes: 1; No: 0.
change_glob_pos = 1;

%% Compute the solution using MPM
[displacement_mpm, velocity_mpm] = MPM_1D(density,...
    Youngs_modulus, gravitational_acceleration, load, height,...
    number_elements, element_size, number_particles_per_element,...
    pos_p_glob, t_step, number_time_steps, total_time, ...
    displacement_p_initial, velocity_p_initial, stress_p_initial,...
    displacement_n_initial, velocity_n_initial, both_ends_fixed, ...
    change_glob_pos, pos_p_loc);

%% Compute the exact solution
w1 = pi*sqrt(Youngs_modulus/density)/height;
b1 = pi/height;

displacement_exact(:,1) = displacement_p_initial;
velocity_exact(:,1) = velocity_p_initial;

for n=1:number_time_steps-1
    for particle = 1:(number_elements*number_particles_per_element)
        velocity_exact(particle,n+1) = 0.1*cos(w1*t_step*n)*...
            sin(b1*pos_p_glob(particle));
        displacement_exact(particle,n+1) = 0.1/w1*sin(w1*t_step*n)*...
            sin(b1*pos_p_glob(particle));
    end
end
clear n particle

%% Flags
% Plot displacement versus time for the selected node? Yes: 1; No: 0 
displ_time = 1; 
% Plot velocity versus time for thr selected node? Yes: 1; No: 0 
velocity_time = 1;

% % Plot displacement versus x-coordinate for the selected node? Yes: 1;
% %No: 0 
% displ_x = 0; 
% % Plot velocity versus x-coordinate for thr selected node? Yes: 1; No: 0 
% velocity_x = 0;
% 
% % Make animation of displacement? Yes: 1; No: 0
% anim_displ = 0;
% 
% % Make animation of velocity? Yes: 1; No: 0
% anim_vel = 0;


%% Plot displacement/velocity versus time for one node
% Select the node
particle_number = floor(number_elements*number_particles_per_element/2); %floor((number_elements*number_particles_per_element)/2);

if displ_time == 1
    figure(1)
    plot_node_displacement_vs_time(particle_number, displacement_mpm...
        (particle_number,:), displacement_exact(particle_number,:),t)
end

if velocity_time == 1
    figure(2)
    plot_node_velocity_vs_time(particle_number,...
        velocity_mpm(particle_number,:),...
        velocity_exact(particle_number,:),t)
end


