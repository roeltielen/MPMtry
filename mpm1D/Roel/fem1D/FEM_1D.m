function [u,v, M_lump, F_int_plot] = FEM_1D(density, E, g, load, H, n_e,...
    element_size, t_step, n_time_steps, u_0, v_0, both_ends_fixed)
% 1D FEM

% INPUT: density, Young's modulus, gravitational acceleration, load, 
% height/length, number of elements, element's size, time step, number of 
% time steps, initial displacement, initial velocity

% OUTPUT: displacement, velocity as function of x and t

%% Set element size
n_n = n_e+1; % number of nodes
h = element_size;


%% Generate element connections 
xvec  = 0:h:H; %global mesh
elements = zeros(n_e,2);
elements_index = zeros(n_e,2);

F_int_plot = zeros(1,n_time_steps);
F_damp = zeros(n_e+1,n_time_steps);

for i = 1:n_e
   elements(i,:) = [xvec(i) xvec(i+1)]; 
   elements_index(i,:) = [i i+1]; %contains the indices of the
    %nodes comprising the elements
end

%% Compute the boolean matrices based on element connections
T = zeros(2, n_n, n_e);
for j=1:n_e
    T(1,elements_index(j,1),j) = 1;
    T(2,elements_index(j,2),j) = 1;
end

%% Initialise the discretized matrices
%..Apply Gaussian rule on the local domain with one Gauss point per element
[M_loc, F_grav_loc, F_int_loc, B_loc] = gauss_rule(density, E, h, g);

%..Assemble the global matrices/vectors created with Gauss rule
[M, F_grav, F_int] = assemble(n_e, T, M_loc, F_grav_loc,...
    F_int_loc);

M_lump = sparse(diag(sum(M,1)));

%..Mass lumping and taking the inverse
M_lump_inv = sparse(diag(1./sum(M,1)));

%..Construct the global traction vector
F_trac = zeros(n_n,1);
F_trac(end) = load;

%% Apply initial conditions
v(:,1) = v_0;
u(:,1) = u_0;
% stress(:,1) = zeros(n_e,1);


%% Time integration
for n=1:n_time_steps-1
    %    stress(:,n) = stress_computation(n_e, E, B_loc, u(:,n));
    F_internal = F_int*u(:,n);
    F_internal;
    
    alpha = 0;
    F_damp = - sign(v(:,n))*alpha.*abs(F_trac + F_grav - F_internal);
    F = F_trac + F_grav - F_internal + F_damp;
    a = M_lump_inv*F;
    
    % Determine internal force at time t = s*dt
    node = floor((n_e+1)/2);
    F_int_plot(1,n) = F_internal(node);
    
    if both_ends_fixed == 1
        a(1) = 0;
        a(n_n) = 0;
    else
        a(1) = 0;
    end
    v(:,n+1) = v(:,n)+ a*t_step;
    if both_ends_fixed == 1
        v(1,n+1) = 0;
        v(n_n,n+1) = 0;
    else
        v(1,n+1) = 0;
    end
    u(:,n+1) = u(:,n) + v(:,n+1)*t_step;
    if both_ends_fixed == 1
        u(1,n+1) = 0;
        u(n_n, n+1) = 0;
    else
        u(1, n+1) = 0;
    end
    
   
  
   clear F_internal F
   
end


