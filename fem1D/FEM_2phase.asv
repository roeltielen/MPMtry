% This code performs 2-phase FEM computation in 1D. It requires an input
% file containing 
% 1) material parameters: porosity, density of the water phase, density of 
% the solid phase, Youngs_modulus, bulk_modulus, Darcy permeability, total
% stress, initial pore pressure, height/length of the sample, gravitational
% acceleration;
% 2) mesh properties: number of elements;
% 3) time descritization properties: total time, CFL number;
% 4) initial conditions: initial velocity of water, initial velocity of
% solid skeleton, initial total stress, initial pore pressure;
% 5) boundary conditions.
% Author: Lisa Wobbes


function [pore_pressure, effective_stress] = FEM_2phase...
    (porosity, density_w,density_s, Youngs_modulus, bulk_modulus,...
    permeability, total_stress, pore_pressure, eff_stress, height,...
    grav_accel, n_e, CFL_number, total_time, velocity_s_initial,...
    velocity_w_initial, total_stress_initial, pore_pressure_initial,...
    eff_stress_initial, v_w_bottom, v_s_bottom, water_traction)


%% Pre-processing

% Compute saturated density and effective stress
density_sat = (1-porosity)*density_s + porosity*density_w;

% Construct mesh
element_size = height/n_e; 
mesh  = 0:element_size:height;

% Waves' velocities
undrained_constr = Youngs_modulus + bulk_modulus/porosity;
vel_undrained = sqrt(undrained_constr/density_sat);
damped_factor = sqrt((porosity*Youngs_modulus/bulk_modulus)/...
    (1-porosity+porosity*Youngs_modulus/bulk_modulus));
vel_damped = damped_factor*sqrt(bulk_modulus/density_w);

% Time step 
t_cr = min(element_size/vel_undrained,element_size/vel_damped);
t_step = 1E-6; %CFL_number*t_cr;
number_time_steps = floor(total_time/t_step); 
t = 0:t_step:(number_time_steps-1)*t_step;

% Generate element connections 
xvec  = 0:element_size:height; %global mesh
elements = zeros(n_e,2);
elements_index = zeros(n_e,2);
for i = 1:n_e
   elements(i,:) = [xvec(i) xvec(i+1)]; 
   elements_index(i,:) = [i i+1]; %contains the indices of the
    %nodes comprising the elements
end

% Compute the boolean matrices based on element connections
n_n = n_e+1; 
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
B_loc = [-1/element_size, 1/element_size];

%..Generate local matrices
M_w_loc = weight*N_loc'*density_w*N_loc*element_size; %m
F_w_grav_loc = [0; 0]; % -weight*N_loc'*density_w*grav_accel*element_size;
Q_loc = weight*N_loc'*(porosity*density_w*grav_accel*element_size...
    /permeability)*N_loc; %m
M_s_loc = weight*N_loc'*(1-porosity)*density_s*N_loc*element_size;%m
M1_w_loc = weight*N_loc'*porosity*density_w*N_loc*element_size;%m
F_grav_loc = [0; 0]; %-weight*N_loc'*density_sat*grav_accel*element_size;

% bottleneck: vereist veel tijd 
%..Assemble matrices/vectors
M_w = assemble_matrix(n_e,T,M_w_loc);
Q = assemble_matrix(n_e,T,Q_loc);
M_s = assemble_matrix(n_e,T,M_s_loc);
M1_w = assemble_matrix(n_e,T,M1_w_loc);
F_w_grav = assemble_vector(n_e,T,F_w_grav_loc);
F_grav = assemble_vector(n_e,T,F_grav_loc);

%..Traction forces
F_w_trac = zeros(n_n,1);
F_w_trac(end) = water_traction; %pore_pressure;
F_trac = zeros(n_n,1);
F_trac(end) = total_stress;

%..Lump the matrices
M_w = lump(M_w);
Q = lump(Q);
M_s = lump(M_s);
M1_w = lump(M1_w);


%..Apply initial conditions
v_w(:,1) = velocity_w_initial;
v_s(:,1) = velocity_s_initial;
pp(:,1) = pore_pressure_initial;
efs(:,1) = eff_stress_initial;


%% Time integration
for nt=1:number_time_steps-1
    nt
%     v_w(1,nt)
%     v_s(1,nt)
    %..Initial internal forces
    F_w_int = zeros(n_n,1);
    for ne = 1:n_e
        F_w_int(elements_index(ne,1)) = F_w_int(elements_index(ne,1)) + ...
            B_loc(1)*pp(ne,nt)*element_size;
        F_w_int(elements_index(ne,2)) = F_w_int(elements_index(ne,2)) + ...
            B_loc(2)*pp(ne,nt)*element_size;
    end
    clear ne

    F_int = zeros(n_n,1);
    for ne = 1:n_e
        F_int(elements_index(ne,1)) = F_int(elements_index(ne,1)) + ...
            B_loc(1)*efs(ne,nt)*element_size +...
            B_loc(1)*pp(ne,nt)*element_size;
        F_int(elements_index(ne,2)) = F_int(elements_index(ne,2)) + ...
            B_loc(2)*efs(ne,nt)*element_size +...
            B_loc(2)*pp(ne,nt)*element_size;
    end
    clear ne

    %..Initial drag force
    %diag(Q)
    F_w_drag = Q*(v_w(:,nt) - v_s(:,nt));
    
    %..Total forces
    F_w = F_w_trac - F_w_int + F_w_grav - F_w_drag;
    F = F_trac - F_int + F_grav;
    
    %..Water acceleration
    a_w = M_w\F_w;
    %a_w(1) = 0;
    
    %..Water velocity
    v_w(:,nt+1) = v_w(:,nt) + a_w*t_step;
    
    %..Solid phase velocity
    v_s(:,nt+1) = v_s(:,nt) + M_s\(-M1_w*a_w + F)*t_step;
    
    %..Boundary conditions
    v_w(1,nt+1) = v_w_bottom;
    v_s(1,nt+1) = v_s_bottom;

    %..Pore pressure increment
    dpp = zeros(n_e,1);
    for ne = 1:n_e
    dpp(ne) = t_step*(bulk_modulus/porosity)*B_loc*((1-porosity)*...
        [v_s(elements_index(ne,1),nt+1); v_s(elements_index(ne,2),nt+1)]... 
        + porosity*[v_w(elements_index(ne,1),nt+1);...
        v_w(elements_index(ne,2),nt+1)]); 
    end
    clear ne
    
    %..Update pore pressure
    pp(:,nt+1) = pp(:,nt) + dpp;
    
    %..Effective stress increment
    defs = zeros(n_e,1);
    for ne = 1:n_e
        defs(ne) = t_step*Youngs_modulus*B_loc*...
            [v_s(elements_index(ne,1),nt+1);...
            v_s(elements_index(ne,2),nt+1)];
    end
    clear ne
    
    %..Update effective stress
    efs(:,nt+1) = efs(:,nt) + defs;
end
clear nt

pore_pressure = pp;
effective_stress = efs;

disp('Yaaaay!')



end