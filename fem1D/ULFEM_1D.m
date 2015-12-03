function [u,v, M_lump, F_int_plot,strain, E_kin, E_pot, E_trac,E_grav] = ULFEM_1D(density, E, g, load, H, n_e,...
    element_size, t_step, n_time_steps, u_0, v_0, both_ends_fixed,mesh,E_kin, E_pot, E_trac,E_grav)
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
strain = zeros(n_e,n_time_steps)

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

%% Apply initial conditions
v(:,1) = v_0;
u(:,1) = u_0;
% stress(:,1) = zeros(n_e,1);
F_damp = zeros(1,n_e+1);


%% Time integration
for n=1:n_time_steps-1
    
    if n > 8
       time = n*t_step; 
    end
    % Initialise the discretized matrices
    %..Apply Gaussian rule on the local domain with one Gauss point per element
    for i = 1:n_e
        h = elements(i,2) - elements(i,1);
        [M_loc(:,:,i), F_grav_loc(:,:,i), F_int_loc(:,:,i), B_loc(:,:,i)] = gauss_rule(density, E, h, g);  
    end
    
    %..Assemble the global matrices/vectors created with Gauss rule
    [M, F_grav, F_int] = assemble_ULFEM(n_e, T, M_loc, F_grav_loc,...
        F_int_loc);

    M_lump = sparse(diag(sum(M,1)));

    %..Mass lumping and taking the inverse
    M_lump_inv = sparse(diag(1./sum(M,1)));

    %..Construct the global traction vector
    F_trac = zeros(n_n,1);
    
    if  n >= 1 && n <= 100 
        F_trac(end) = (n/100)*load;
    else 
        F_trac(end) = load;
    end
        
    
    %    stress(:,n) = stress_computation(n_e, E, B_loc, u(:,n));
    F_internal = F_int*u(:,n);
    F_internal;
    
    alpha = 0;
    F_damp = - sign(v(:,n))*alpha.*abs(F_trac+F_grav-F_internal);
    F = F_trac + F_grav - F_internal + F_damp;
    a = M_lump_inv*F;
    
    % Determine internal force at time t = s*dt
    node = floor((n_e+1)/2);
    %node = 2;
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
   
   
    
   % Update position of the nodes
   elements(:,1) = elements(:,1) + v(1:n_e,n+1)*t_step;
   elements(:,2) = elements(:,2) + v(2:n_e+1,n+1)*t_step;
    
   % FOR LOOP
   % Determine B  for every h 
   % epsilon_i(n+1) = epsilon_i(n) + v(n+1)*t_step*B(n)
   % END 
   % Strain in Gauss Points
   for i = 1:n_e 
      h = elements(i,2) - elements(i,1); 
      B_temp = [-1/h , 1/h]; 
      strain(i,n+1) = strain(i,n) + B_temp*[v(i,n+1)*t_step ; v(i+1,n+1)*t_step];
   end 
   
    
    % Determine the kinetic and potential energy of the system
    for s = 1:n_n
       E_kin(1,n+1) = E_kin(1,n+1) + 0.5*M_lump(s,s)*v(s,n+1)*v(s,n+1); 
       E_pot(1,n+1) = E_pot(1,n+1) + -0.5*F_internal(s)*u(s,n+1);
       E_trac(1,n+1) = E_trac(1,n+1) + F_trac(s)*u(s,n+1);
       E_grav(1,n+1) = E_grav(1,n+1) + F_grav(s)*u(s,n+1);
    end
    
   clear F_internal F
end

