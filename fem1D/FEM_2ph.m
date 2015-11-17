% Oedometer 1D

% This file provides input for FEM to computes numerical solution, provides
% an exact solution, illustrates results and provides the global error

clear all
close all
beep off
clc 

%% Input
% Constants
porosity = 0.4;
density_w = 1E3;
density_s = 2.6E3;
density_sat = (1-porosity)*density_s + n*density_w;
Youngs_modulus = 5E9; %check
bulk_modulus = 2E9; %check
permeability = 1E-5; %check
grav_accel = -9.8; 
total_stress = 1E3; %check
pore_pressure = 1E3; %check
height = 2.5; %height/length

% Mesh properties
number_elements = 16; % number of elements
element_size = height/number_elements; 
mesh  = 0:element_size:height; %mesh: NEEDED FOR INITIAL VELOCITY AND
%EXACT SOLUTION

% Waves' velocities
undrained_constr = Youngs_modulus + bulk_modulus/porosity;
vel_undrained = sqrt(undrained_constr/density_sat);
damped_factor = sqrt((porosity*Youngs_modulus/bulk_modulus)/...
    (1-porosity+porosity*Youngs_modulus/bulk_modulus));
vel_damped = damped_factor*sqrt(bulk_modulus/density_w);

% Time step %check!
CFL_number = 0.9;
total_time = 1; 
t_cr = min(h/vel_undrained,h/vel_damped);
t_step = CFL_number*t_cr;
number_time_steps = floor(total_time/t_step); % set here the total time
t = 0:t_step:(number_time_steps-1)*t_step;

% Initial conditions 
velocity_s_initial = zeros(number_elements + 1,1);
velocity_w_initial = zeros(number_elements + 1,1);
total_stress_vec = total_stress*ones(number_elements + 1,1);
pore_pressure_vec = pore_pressure*ones(number_elements + 1,1);

% Boundary conditions
both_ends_fixed = 0;

%% Compute the solution using FEM
% Set element size
n_n = n_e+1; % number of nodes
h = element_size;

% Generate element connections 
xvec  = 0:h:H; %global mesh
elements = zeros(n_e,2);
elements_index = zeros(n_e,2);
for i = 1:n_e
   elements(i,:) = [xvec(i) xvec(i+1)]; 
   elements_index(i,:) = [i i+1]; %contains the indices of the
    %nodes comprising the elements
end

% Compute the boolean matrices based on element connections
T = zeros(2, n_n, n_e);
for j=1:n_e
    T(1,elements_index(j,1),j) = 1;
    T(2,elements_index(j,2),j) = 1;
end

% Gaussian rule on the local domain with one Gauss point per element
pos_loc = 0.5;
weight  = 1;

%..Define basis functions vector 
N_loc = [1-pos_loc, pos_loc];

%..Define strain-displacement matrix 
B_loc = [-1/h, 1/h];

%..Local matrices
M_w_loc = weight*N_loc'*density_w*N_loc*h;
F_w_grav_loc = -weight*N_loc*density_w*grav_accel*h;
Q_loc = weight*N_loc'*porosity*density_w*grav_accel*h/permeability;
M_s_loc = weight*N_loc'*(1-porosity)*density_s*N_loc*h;
M1_w_loc = weight*N_loc'*porosity*density_w*N_loc*h;
F_grav_loc = -weight*N_loc'density_sat*grav_accel*h;

