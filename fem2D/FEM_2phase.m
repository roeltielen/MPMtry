% This code performs 2-phase FEM computation in 2D for isotropic materials.
% It requires an input file containing: 
% 1) material parameters: porosity, density of the water phase, density of 
% the solid phase, Youngs_modulus, bulk_modulus, Darcy permeability, total
% stress, initial pore pressure, height/length of the sample, width of...
% the sample gravitational acceleration;
% 2) mesh properties: number of elements in x direction; number of elements
% in y direction
% 3) time descritization properties: total time, CFL number;
% 4) initial conditions: initial velocity of water, initial velocity of
% solid skeleton, initial total stress, initial pore pressure;
% 5) boundary conditions.
% Author: Lisa Wobbes

function [] = FEM_2phase...
    (porosity, density_w,density_s, Youngs_modulus, bulk_modulus,...
    permeability, total_stress, pore_pressure, eff_stress, height, width,...
    grav_accel, n_ey, n_ex, CFL_number, total_time, velocity_s_initial,...
    velocity_w_initial, total_stress_initial, pore_pressure_initial,...
    eff_stress_initial, v_w_bottom, v_s_bottom, water_traction)
%% Pre-processing

% Compute saturated density and effective stress
density_sat = (1-porosity)*density_s + porosity*density_w;

% Construct mesh
element_length = height/n_ey; 
element_width = width/n_ex; 
yvec  = [0:element_length:height]
xvec  = [0:element_width:width]
% X and Y contain x- and y- position of each node
[X,Y] = meshgrid(xvec,yvec)

% Waves' velocities
undrained_constr = Youngs_modulus + bulk_modulus/porosity;
vel_undrained = sqrt(undrained_constr/density_sat);
damped_factor = sqrt((porosity*Youngs_modulus/bulk_modulus)/...
    (1-porosity+porosity*Youngs_modulus/bulk_modulus));
vel_damped = damped_factor*sqrt(bulk_modulus/density_w);

% Time step 
t_cr = min(element_length/vel_undrained,element_length/vel_damped);
t_step = 1E-6; %CFL_number*t_cr;
number_time_steps = floor(total_time/t_step); 
t = 0:t_step:(number_time_steps-1)*t_step;

% % Generate element connections 
% elements = zeros(n_e,2);
% elements_index = zeros(n_e,2);
% for i = 1:n_e
%    elements(i,:) = [xvec(i) xvec(i+1)]; 
%    elements_index(i,:) = [i i+1]; %contains the indices of the
%     %nodes comprising the elements
% end

end