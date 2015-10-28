% Oedometer with load, without self-weight

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
load = -5E3;
height = 1; %height/length

% Mesh properties
number_elements = 8; % number of elements
element_size = height/number_elements; 
mesh  = 0:element_size:height; %mesh: NEEDED FOR INITIAL VELOCITY AND
%EXACT SOLUTION

% Time step 
format long
CFL_number = 1;
total_time = 1.5; 
t_cr = element_size/sqrt(Youngs_modulus/density);
t_step = CFL_number*t_cr;
number_time_steps = floor(total_time/t_step); % set here the total time
t = 0:t_step:(number_time_steps-1)*t_step;




%% Compute an exact solution

for n=1:number_time_steps
    for node = 1:number_elements+1
        [position_exact(node,n),displacement_exact(node,n),...
            velocity_exact(node,n)] = exact_solution...
            (density,Youngs_modulus, load,-gravitational_acceleration,...
            height,mesh(node), t_step*(n-1));
    end
end
clear n node

% T_step = 51;
% T = T_step*t_step
% T_step = T_step+1; 

T = 1;
T_step = floor(T/t_step);

for node = 1:number_elements + 1
    [position_exactT(node,:),displacement_exactT(node,:),...
        velocity_exactT(node,:)] = exact_solution...
        (density,Youngs_modulus, load,-gravitational_acceleration,...
        height,mesh(node), T);
end
clear node
displacement_exactT

%% % Initial conditions 
displacement_initial = displacement_exact(:,1);
velocity_initial = velocity_exact(:,1);

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
% T = 0.1;
% T_step = floor(T/t_step);

if displ_x == 1
    figure(3)
    plot_displacement_vs_x(T, displacement_fem(:,T_step),...
        displacement_exactT, mesh)
end

if velocity_x == 1
    figure(4)
    plot_velocity_vs_x(T, velocity_fem(:,T_step),...
        velocity_exactT, mesh)
end


%% Animation displacement/velocity versus x-coordinate
time_interval = 1:4:number_time_steps;

if anim_displ == 1
    for i = 1:length(time_interval)
        n = time_interval(i);

        max_displ = max(displacement_fem(floor((number_elements+1)/2),:));
        max_displ = max_displ*1.1;

        min_displ = min(displacement_fem(floor((number_elements+1)/2),:));
        min_displ = min_displ*1.1;

        figure(5)
        animation_displacement_vs_x(n*t_step, displacement_fem(:,n),...
            displacement_exact(:,n), mesh, max_displ, min_displ, height)
        pause(0.1);
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
displacement_fem(:,T_step)



