% 1D FEM for Oedometer problem
clear all
close all
clc

%..Constants
density = 1E3;
E = 1E5; % Young's modulus
g = -9.8; % gravitational acceleration
load = 0; 

%..Set number of elements and number of nodes (inclusing the boundary nodes)
n_e = 4;
n_n = n_e+1; 

%..Generate mesh size and vector of grid points.. 
H = 1; %height of the column
h = H/n_e; % element size
position = 1; %required for the plotting

%..Generate global mesh, element connections 
xvec  = 0:h:1; %global mesh
elements = zeros(n_e,2);
for i = 1:n_e
   elements(i,:) = [xvec(i) xvec(i+1)]; 
end

%..Compute the boolean matrices based on element connections
T = zeros(2, n_n, n_e);
for j=1:n_e
   node_1 = elements(j,1);
   node_2 = elements(j,2);
   pos_1 = round(node_1/h)+1;
   pos_2 = round(node_2/h)+1;
   T(1,pos_1,j) = 1;
   T(2,pos_2,j) = 1;
   clear node_1 node_2 pos_1 pos_2
end

%..Time step 
CN = 0.5;
t_cr = h/sqrt(E/density);
t_step = CN*t_cr
n_time_steps = 2.5/t_step; % whatever
t = 0:t_step:(n_time_steps-1)*t_step;

%..Apply Gaussian rule on the local domain with one Gauss point per element
[M_loc, F_grav_loc, F_int_loc, B_loc] = gauss_rule(density, E, h, g);

%..Assemble the global matrices/vectors created with Gauss rule
[M, F_grav, F_int] = assemble (n_e, T, M_loc, F_grav_loc,...
    F_int_loc);

%..Mass lumping and taking the inverse
M_lump_inv = sparse(diag(1./sum(M,1)));

%..Construct the global traction vector
F_trac = zeros(n_n,1);
F_trac(end) = load;

%..Initial conditions
v(:,1) = zeros(n_n,1);
u(:,1) = zeros(n_n,1);
stress(:,1) = zeros(n_e,1);

%..Time integration
for n=1:n_time_steps-1
  % stress(:,n) = stress_computation(n_e, E, B_loc, u(:,n));
   F_internal = F_int*u(:,n);
   F = F_trac + F_grav - F_internal;
   v(:,n+1) = v(:,n)+ M_lump_inv*F*t_step;
   u(:,n+1) = u(:,n) + v(:,n+1)*t_step;
   u(1,n+1) = 0; %boundary condition
   clear F_internal F
end


%.. Exact solution
exact_sol = exact_solution_original_oedometer(density, E, load, -g, H, position, t);

%..Plot
figure;
plot(t, position+u(end,:))
hold on
plot(t', exact_sol, 'r')

%..Plot stress
xgauss = xvec+0.5*h
xgauss(end) = [];
size(stress(:,4))
size(xgauss)
figure;
plot(stress(:,round(0.15/t_step)), xgauss);
set(gca, 'Xdir', 'reverse')



