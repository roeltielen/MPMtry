% OLD: one-dimensional wave propagation
clear all
close all
beep off
clc 

%% Input
% Constants
porosity = 0.4;
density_w = 1E3;
density_s = 2.6E3;
density_sat = (1-porosity)*density_s + porosity*density_w;
Youngs_modulus = 2.5E9; %check
bulk_modulus = 2E9; %check
permeability = 1E-5; %check
grav_accel = 10; %10; 
total_stress = -10E3; %check
pore_pressure = -10E3; %check
height = 2; %height/length
eff_stress = total_stress - pore_pressure;

% Mesh properties
n_e = 800; % number of elements
element_size = height/n_e; 
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
total_time = 3E-3; %0.0018; 
t_cr = min(element_size/vel_undrained,element_size/vel_damped);
t_step = 1E-6; %CFL_number*t_cr;
number_time_steps = floor(total_time/t_step); 
t = 0:t_step:(number_time_steps-1)*t_step;

% Initial conditions 
velocity_s_initial = zeros(n_e + 1,1);
velocity_w_initial = zeros(n_e + 1,1);
total_stress_initial = zeros(n_e,1);
total_stress_initial(end) = total_stress;
pore_pressure_initial = zeros(n_e,1);
%pore_pressure_initial(end) = pore_pressure;
eff_stress_initial = zeros(n_e,1);
%eff_stress_initial(end) = eff_stress;

%% Compute the solution using FEM
% Set element size
n_n = n_e+1; % number of nodes

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

%..Assemble matrices/vectors
M_w = assemble_matrix(n_e,T,M_w_loc);
Q = assemble_matrix(n_e,T,Q_loc);
M_s = assemble_matrix(n_e,T,M_s_loc);
M1_w = assemble_matrix(n_e,T,M1_w_loc);
F_w_grav = assemble_vector(n_e,T,F_w_grav_loc);
F_grav = assemble_vector(n_e,T,F_grav_loc);

%..Traction forces
F_w_trac = zeros(n_n,1);
F_w_trac(end) = pore_pressure;
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

%..Time integration
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
    v_w(1,nt+1) = 0;
    v_s(1,nt+1) = 0;

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

% analytical solution for stress
[t1, AS, AP] = as_wave_propagation(permeability);

figure(2);
plot(t,efs(floor(3*n_e/4),:)/total_stress,'LineWidth',2)
elements(floor(3*n_e/4),:)
hold on
plot(t1,AS,'--k','LineWidth',2)
xlabel('time [s]', 'FontSize', 12)
set(gca,'FontSize',11)
ylabel('normalized effective stress [-]','FontSize', 12)
%title(sprintf('Displacement at time %g [s]',T)) 
legend('FEM', 'Exact', 'Location','southeast')
% hold on
set(0,'DefaultFigureColor',[1 1 1])
set(gcf, 'PaperPosition', [0 0 6 6]);
set(gcf, 'PaperSize', [6 6]);


figure(3);
plot(t,pp(floor(3*n_e/4),:)/total_stress,'r','LineWidth',2)
hold on
plot(t1,AP,'--k','LineWidth',2)
xlabel('time [s]', 'FontSize', 12)
set(gca,'FontSize',11)
ylabel('normalized pore presssure [-]','FontSize', 12)
% title(sprintf('Displacement at time %g [s]',T)) 
legend('FEM', 'Exact', 'Location','southeast')
% hold on
set(0,'DefaultFigureColor',[1 1 1])
set(gcf, 'PaperPosition', [0 0 6 6]);
set(gcf, 'PaperSize', [6 6]);

% figure(4);
% x = element_size/2:element_size:n_e*element_size;
% for nt = 1:20:number_time_steps-1
%    plot(x',pp(:,nt)) 
%    pause(1)
% end






