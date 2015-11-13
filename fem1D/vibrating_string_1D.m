% Vibrating string with fixed ends

% This file provides input for FEM to computes numerical solution, computes
% an exact solution, illustrates results and provides the global error

clear all
close all
beep off
clc 

%% Input
% Constants
density = 1E3; %1;
Youngs_modulus = 1E5; %100;
gravitational_acceleration = 0; 
load = 0;
height = 1; %height/length

% Mesh properties
number_elements = 4; % number of elements
element_size = height/number_elements; 
mesh  = 0:element_size:height; %mesh: NEEDED FOR INITIAL VELOCITY AND
%EXACT SOLUTION

% Time step 
CFL_number = 0.9;
total_time = 1.5; 
t_cr = element_size/sqrt(Youngs_modulus/density);
t_step = 1E-4; %CFL_number*t_cr;
number_time_steps = floor(total_time/t_step); % set here the total time
t = 0:t_step:(number_time_steps-1)*t_step;

% Initial conditions 
displacement_initial = zeros(number_elements + 1,1);
velocity_initial = 0.1*sin(pi/height*mesh);

% Boundary conditions
both_ends_fixed = 1;

%% Compute the solution using FEM
[displacement_fem,velocity_fem, M_lump] = FEM_gt(density,...
    Youngs_modulus, gravitational_acceleration, load, height,...
    number_elements, element_size, t_step, number_time_steps,...
    displacement_initial, velocity_initial, both_ends_fixed);

%% Compute the exact solution
w1 = pi*sqrt(Youngs_modulus/density)/height;
b1 = pi/height;

displacement_exact(:,1) = displacement_initial;
velocity_exact(:,1) = velocity_initial;

for n=1:number_time_steps-1
    for node = 1:number_elements + 1
        velocity_exact(node,n+1) = 0.1*cos(w1*t_step*n)*sin(b1*mesh(node));
        displacement_exact(node,n+1) = 0.1/w1*sin(w1*t_step*n)*...
            sin(b1*mesh(node));       
    end
end
clear n node

% Solution at certain time (required for accuracy)
% T = 0.0003;
% T_step = floor(T/t_step)+1
T_step = 510;
T = T_step*t_step;
T_step = T_step+1;

for node = 1:number_elements + 1
    velocity_exactT(node,:) = 0.1*cos(w1*T)*sin(b1*mesh(node));
    displacement_exactT(node,:) = 0.1/w1*sin(w1*T)*sin(b1*mesh(node));
end
clear node
%displacement_exactT

%% Flags
% Plot displacement versus time for the selected node? Yes: 1; No: 0 
displ_time = 0; 
% Plot velocity versus time for thr selected node? Yes: 1; No: 0 
velocity_time = 0;

% Plot displacement versus x-coordinate for the selected node? Yes: 1;
%No: 0 
displ_x = 0; 
% Plot velocity versus x-coordinate for thr selected node? Yes: 1; No: 0 
velocity_x = 0;

% Make animation of displacement? Yes: 1; No: 0
anim_displ = 0;

% Make animation of velocity? Yes: 1; No: 0
anim_vel = 0;

% Compute the global error at time T_step? Yes: 1, No: 0
compute_error = 1;


%% Plot displacement/velocity versus time for one node
% Select the node
node_number = floor((number_elements+1)/2);

if displ_time == 1
    figure(1)
    plot_node_displacement_vs_time(node_number, displacement_fem...
        (node_number,:), displacement_exact(node_number,:),t)
end

if velocity_time == 1
    figure(2)
    plot_node_velocity_vs_time(node_number, velocity_fem(node_number,:),...
        velocity_exact(node_number,:),t)
end

%% Plot displacement/velocity versus x-coordinate for a certain moment
if displ_x == 1
    figure(3)
    plot_displacement_vs_x(T, displacement_fem(:,T_step),...
        displacement_exactT, mesh)
end

if velocity_x == 1
    figure(4)
    plot_velocity_vs_x(T, velocity_fem(:,T_step), velocity_exactT, mesh)
end


%% Animation displacement/velocity versus x-coordinate
time_interval = 1:150:15000;

if anim_displ == 1
    for i = 1:length(time_interval)
        n = time_interval(i);

        max_displ = max(displacement_fem(floor((number_elements+1)/2),:));
        max_displ = max_displ*1.1;

        min_displ = min(displacement_fem(floor((number_elements+1)/2),:));
        min_displ = min_displ*1.1;

        figure(5)
        animation_displacement_vs_x(n*t_step, displacement_fem(:,n),...
            displacement_exact(:,n), mesh, 0.004, -0.004, height)
        pause(0.3);
    end
end


if anim_vel == 1
    for i = 1:length(time_interval)
        n = time_interval(i);

        max_vel = max(velocity_fem(floor((number_elements+1)/2),:));
        max_vel = max_vel*1.1;

        min_vel = min(velocity_fem(floor((number_elements+1)/2),:));
        min_vel = min_vel*1.1;

        figure(6)
        animation_velocity_vs_x(n*t_step, velocity_fem(:,n),...
            velocity_exact(:,n), mesh, max_vel, min_vel)
        pause(0.1);
    end
end

%% Error computation

if compute_error == 1
    % Compute the norm of the discretization error in 2-norm

    errnrm = compute_error_norm(displacement_fem(:,T_step),...
        displacement_exactT, M_lump);
    fprintf('For time t = %e\n', T)
    fprintf('The error norm = %e\n', errnrm)
end

%  displacement_exactT
