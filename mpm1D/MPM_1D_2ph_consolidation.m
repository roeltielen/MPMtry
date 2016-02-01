% Oedometer 1D

% This file provides input for FEM to computes numerical solution, provides
% an exact solution, illustrates results and provides the global error

clear all
close all
beep off
clc

tic
%% Input
% Constants
porosity = 0.4;
density_w = 1E3;
density_s =  2600; %(1/0.6)*1.6E3;
density_sat = (1-porosity)*density_s + porosity*density_w;
Youngs_modulus = 2.5E9; %1E7; %check
bulk_modulus = 2E9; %3E8; %check
permeability = 1E-5; %check %hydraulic conductivity is a better name
hydraulic_conductivity = 1E-5;
grav_accel = 10;
total_stress = -10E3; %check
pore_pressure = -10E3; %check
height = 1; %height/length
eff_stress = total_stress - pore_pressure;

cv = permeability/(density_w*grav_accel*(1/Youngs_modulus+...
    porosity/bulk_modulus));
Tf = [0.02]; %, 0.05, 0.1, 0.2, 0.5, 1.0];

% Mesh properties
n_e = 300; % number of elements
element_size = height/n_e;
mesh  = 0:element_size:height;

% Waves' velocities
undrained_constr = Youngs_modulus + bulk_modulus/porosity;
vel_undrained = sqrt(undrained_constr/density_sat);
damped_factor = sqrt((porosity*Youngs_modulus/bulk_modulus)/...
    (1-porosity+porosity*Youngs_modulus/bulk_modulus));
vel_damped = damped_factor*sqrt(bulk_modulus/density_w);

% Time step %check!
CFL_number = 0.9;
total_time = 0.012005; %0.600005;
t_cr = min(element_size/vel_undrained,element_size/vel_damped);
t_step = 1E-6; %CFL_number*t_cr;
n_time_steps = floor(total_time/t_step);
t = 0:t_step:(n_time_steps-1)*t_step;


%%
%%%Compute the solution using MPM
%% Mesh and element connections
% Set element size
n_n = n_e+1; % number of nodes

% Generate element connections
xvec  = 0:element_size:height; %global mesh
elements = zeros(n_e,2);
elements = sparse(elements);
elements_index = zeros(n_e,2);
elements = sparse(elements);
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

%% Particle initialization
% Number of particles
n_ep = 2;
n_p = n_e*n_ep;

% Initial particle location within elements
elements_particles = zeros(n_e,n_p);

for i=1:n_e
    for j = (i-1)*n_ep+1:i*n_ep
        elements_particles(i,j) = 1;
    end
end
clear i

% Particles are equally distributed over element.
% initial local position of particles
pos_p_loc = ones(n_e*n_ep,1)/(n_ep+1);
for i = 2:n_ep
    pos_p_loc(i:n_ep:n_e*n_ep,1) = i/(n_ep+1);
end

% initial global position of particles
pos_p_minloc = zeros(n_e*n_ep,1);
for i = 1:n_e-1
    pos_p_minloc(i*n_ep+1:(i+1)*n_ep,1)...
        = pos_p_minloc(i*n_ep+1:(i+1)*n_ep,1) + i*element_size;
end
pos_p_glob = pos_p_minloc + element_size*pos_p_loc;

% Initial particle properties
volume_p = (1/n_ep)*element_size*ones(n_p,1);
density_w_p = density_w*ones(n_p,1);
density_s_p = density_s*ones(n_p,1);
density_sat_p = density_sat*ones(n_p,1);
porosity_p = porosity*ones(n_p,1);
hydr_conduct_p = hydraulic_conductivity*ones(n_p,1);
bulk_mod_p = bulk_modulus*ones(n_p,1);
m_s_p = density_s_p.*volume_p;
m_w_p = density_w_p.*volume_p;
m_p = density_sat_p.*volume_p;
f_grav_s_p = m_s_p*grav_accel;
f_grav_w_p = m_w_p*grav_accel;
f_grav_p = m_p*grav_accel;
f_trac_p = zeros(n_p,1);
f_trac_w_p = zeros(n_p,1);
%f_trac_w_p(end) = pore_pressure;

% Initial conditions
velocity_s_initial = zeros(n_p,1);
velocity_w_initial = zeros(n_p,1);
total_stress_initial = total_stress*ones(n_p,1);
pore_pressure_initial = pore_pressure*ones(n_p,1);
eff_stress_initial = zeros(n_p,1);

% Vectors to store information
v_s_p = zeros(n_p, n_time_steps-1);
v_s_p(:,1) = velocity_s_initial;
v_w_p = zeros(n_p, n_time_steps-1);
v_w_p(:,1) = velocity_w_initial;
pore_pressure_p = zeros(n_p, n_time_steps-1);
pore_pressure_p(:,1) = pore_pressure_initial;
efs_p = zeros(n_p, n_time_steps-1);
efs_p(:,1) = eff_stress_initial;

%% Time integration
for nt=1:1 %number_time_steps-1
    % Basis functions and derivatives in material points
    N_vec = [ones(n_p,1) - pos_p_loc pos_p_loc];
    B_vec = [-1/element_size*ones(n_p,1) 1/element_size*ones(n_p,1)];

    % Specific boundary conditions for consolidation problem
    f_trac_p(end) = total_stress;

    % Determine the number of particles in each element
    num_particles_in_e = sum(elements_particles,2);

    % Initialisation of element matrix and vectors
    M_w_e = zeros(2,2,n_e);
    M1_w_e = zeros(2,2,n_e);
    M_s_e = zeros(2,2,n_e);
    P_s_e = zeros(2,n_e);
    P_w_e = zeros(2,n_e);
    F_grav_e = zeros(2,n_e); 
    F_grav_w_e = zeros(2,n_e);
    F_int_e = zeros(2,n_e);
    F_int_w_e = zeros(2,n_e);
    Q_e = zeros(2,2,n_e);
    F_trac_w = zeros(2,n_e);
    F_trac = zeros(2,n_e);

    % Determine element (lumped) matrix and vectors for active elements
    for el = 1:n_e
        if num_particles_in_e(el) > 0
            for p = element_particle(elements_particles(el,:))
                % Mass matrices
                M_w_e(1,1,el)  = M_w_e(1,1,el) + m_w_p(p)*N_vec(p,1);
                M_w_e(2,2,el)  = M_w_e(2,2,el) + m_w_p(p)*N_vec(p,2);
                M1_w_e(1,1,el) = M1_w_e(1,1,el) +...
                    porosity_p(p)*m_w_p(p)*N_vec(p,1);
                M1_w_e(2,2,el) = M1_w_e(2,2,el) +...
                    porosity_p(p)*m_w_p(p)*N_vec(p,2);
                M_s_e(1,1,el)  = M_s_e(1,1,el) +...
                   (1 - porosity_p(p))*m_s_p(p)*N_vec(p,1);
                M_s_e(2,2,el)  = M_s_e(2,2,el) +...
                   (1 - porosity_p(p))*m_s_p(p)*N_vec(p,1);
               % Momentum vectors
                P_s_e(:,el) = P_s_e(:,el) +...
                    (1-porosity_p(p))*m_s_p(p)*(N_vec(p,:))'*v_s_p(p,nt); 
                P_w_e(:,el) = P_w_e(:,el) +...
                    porosity_p(p)*m_w_p(p)*(N_vec(p,:))'*v_w_p(p,nt);
               % Force vectors
               F_grav_e(:,el) = F_grav_e(:,el) + (N_vec(p,:))'*f_grav_p(p);
               F_grav_w_e(:,el) = F_grav_w_e(:,el) +...
                   (N_vec(p,:))'*f_grav_p(p);
               F_int_w_e(:,el) = F_int_w_e(:,el) +...
                   (B_vec(p,:))'*pore_pressure_p(p,nt)*volume_p(p);
               % dit moet ook met total stress kunnen!!!!!!!!!!!!!!!!!!
               F_int_e(:,el) = F_int_e(:,el) +...
                   (B_vec(p,:))'*pore_pressure_p(p,nt)*volume_p(p) +...
                   (B_vec(p,:))'*efs_p(p,nt)*volume_p(p);
               % Drag force matrix Q
               % I assumed constant permeability for all particles!!!!!!!!
               Q_e(1,1,el) = Q_e(1,1,el) +... 
                   porosity_p(p)*m_w_p(p)*grav_accel*N_vec(p,1)/...
                   permeability;
               Q_e(2,2,el) = Q_e(2,2,el) +... 
                   porosity_p(p)*m_w_p(p)*grav_accel*N_vec(p,2)/...
                   permeability;
            end
            clear p
        end
    end
    clear el
    
    % Traction forces
    % Traction force for water phase stays unchanged for the consolidation
    % problem
    el = particle_element(elements_particles(:,end)); 
    F_trac_e(:,el) = (N_vec(end,:))'*f_trac_p(end);
    clear el
    
    
    
    
    

end








