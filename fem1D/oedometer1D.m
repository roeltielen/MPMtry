% Oedometer 1D

% This file provides input for FEM to computes numerical solution, provides
% an exact solution, illustrates results and provides the global error

clear all
close all
beep off
clc 

%% Input
% Constants
density = 1E3;
Youngs_modulus = 1E5;
gravitational_acceleration = -9.8; 
load = -5E3;
height = 1; %height/length

% Mesh properties
number_elements = 4; % number of elements
element_size = height/number_elements; 
mesh  = 0:element_size:height; %mesh: NEEDED FOR INITIAL VELOCITY AND
%EXACT SOLUTION

% Time step 
CFL_number = 0.1;
total_time = 1.5; 
t_cr = element_size/sqrt(Youngs_modulus/density);
t_step = 1E-4; %CFL_number*t_cr;
number_time_steps = floor(total_time/t_step); % set here the total time
t = 0:t_step:(number_time_steps-1)*t_step;

% Initial conditions 
displacement_initial = zeros(number_elements + 1,1);
velocity_initial = zeros(number_elements + 1,1);

% Boundary conditions
both_ends_fixed = 0;

%% Compute the solution using FEM
[displacement_fem,velocity_fem, M_lump] = FEM_1D(density,...
    Youngs_modulus, gravitational_acceleration, load, height,...
    number_elements, element_size, t_step, number_time_steps,...
    displacement_initial, velocity_initial, both_ends_fixed);

% Compute the position of nodes by adding the initial position to the
% displacement

position_fem = zeros(size(displacement_fem));

for n = 1:number_time_steps
    position_fem(:,n) = displacement_fem(:,n) + mesh';
end
clear n


%% Obtain the exact solution
position_exact = zeros(number_elements, number_time_steps);

for node = 1:number_elements + 1
    [position_exact(node,:),dummy,dummy] = exact_solution(density,...
        Youngs_modulus, load, -gravitational_acceleration, height,...
        mesh(node), t);
    clear dummy
end
clear node

T = 0.1;
T_step = floor(T/t_step);
for node = 1:number_elements + 1
    [position_exactT(node,:),displacement_exactT(node,:),...
        velocity_exactT(node,:)] = exact_solution...
        (density,Youngs_modulus, load,-gravitational_acceleration,...
        height,mesh(node), T);
end
clear node
displacement_exactT

%% Flags
% Plot displacement versus time for the selected node? Yes: 1; No: 0 
displ_time = 1; 

% Plot displacement versus x-coordinate for the selected node? Yes: 1;
%No: 0 
displ_x = 0; 

% Make animation of displacement? Yes: 1; No: 0
anim_displ = 0;

% Compute the global error at time T_step? Yes: 1, No: 0
compute_error = 1;


%% Plot displacement versus time for one node
% Select the node
node_number = floor((number_elements+1)/2);

if displ_time == 1
    figure(1)
    plot_node_displacement_vs_time(node_number, position_fem...
        (node_number,:), position_exact(node_number,:),t)
end


%% Plot displacement versus x-coordinate for a certain moment
% Select the time moment
% T_step = floor(number_time_steps/10);
% T = t_step*T_step;

T = 0.1;
T_step = floor(T/t_step);

if displ_x == 1
    figure(3)
    plot_displacement_vs_x(T, position_fem(:,T_step),...
        position_exact(:,T_step), mesh)
end


%% Animation displacement versus x-coordinate
time_interval = 1:20:number_time_steps;

if anim_displ == 1
    for i = 1:length(time_interval)
        n = time_interval(i);

        max_displ = 1;

        min_displ = 0;

        figure(5)
        animation_displacement_vs_x(n*t_step, position_fem(:,n),...
            position_exact(:,n), mesh', max_displ, min_displ, height)
        pause(0.1);
    end
end


%% Error computation

if compute_error == 1
    % Compute the norm of the discretization error in 2-norm
    format long
    position_fem(:,T_step)
    position_exact(:,T_step)
    errnrm = compute_error_norm(position_fem(:,T_step),...
        position_exact(:,T_step), M_lump);
    fprintf('For time t = %e\n', T)
    fprintf('The error norm = %e\n', errnrm)
end




