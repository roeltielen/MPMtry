function[u_p, v_p, mass_p,F_int_plot,stress_particle,indicator,v_n] = MPM_1D_Linear_vibrating_string(density, E, g, load,H, l, n_e, h, n_ep,...
    pos_p_glob,  pos_p_loc, t_step, n_time_steps, total_time, u_p, v_p, stress_p, u_n,...
    v_n, both_ends_fixed, change_glob_pos,tau,langranian,reset) 
 % MPM_1D _Linear_vibrating_string(uses element_particle.m, particle_element.m, assemble_matrix,
 % assemble_vector)
 % INPUT: density, Young's modulus, gravitational acceleration, load,
 % height, number of elements, element size, number of particles per
 % element, global particle position, local particle positions, time step, 
 % number of time steps, total time, initial particle displacement, initial 
 % particle velocity, initial particle stress, initial nodal displacement, 
 % initial nodal velocity, bounday conditions, update global positions)
 % OUTPUT: displacement vector, velocity vector and mass vector for each time step



%% Mesh and element connections
% Define nodes and mesh
n_n = n_e+1; 
mesh  = 0:H/n_e:H; 

% Define elements and index of elements
elements_nodes = zeros(n_e,2);
elements_nodes_index = zeros(n_e,2);
for el = 1:n_e
    elements_nodes(el,:) = [mesh(el) mesh(el+1)]; 
    elements_nodes_original(el,:) = [mesh(el) mesh(el+1)]; 
    elements_nodes_index(el,:) = [el el+1]; 
end
clear el

%% Particle initialization
% Number of particles
n_p = n_e*n_ep; 

% Initial particle location
%elements_particles(:,:,1) = zeros(n_e,n_p);

%for i=1:n_e
%   for j = (i-1)*n_ep+1:i*n_ep
%     elements_particles(i,j) = 1;                                     
%   end
%end

% Initial particle location
elements_particles = zeros(n_e,n_p);

 for p = 1:n_p
    el = (pos_p_glob(p) - mod(pos_p_glob(p),h))/h; 
    elements_particles(floor(el+1),p) = 1;
 end
clear p


% Initial particle properties
volume_p = (1/n_ep)*h*ones(n_p,1); 
density_p = density*ones(n_p,1); 
mass_p = volume_p.*density_p; 
f_grav_p = mass_p*g; 
f_trac_p = zeros(n_p,1); 

% History dependent deformation gradient
F_p = zeros(n_p,n_time_steps-1);

v_n = zeros(n_n,n_time_steps);
F_int_plot = zeros(1,n_time_steps);
F_ext_plot = zeros(1,n_time_steps);
F_diff = zeros(1,n_time_steps);
stress_particle = zeros(1,n_time_steps);
indicator = zeros(n_n,n_time_steps);

%% Time integration
for s = 1:n_time_steps-1
    
    % Basis functions and derivatives in material points
    N_vec = [ones(n_p,1) - pos_p_loc pos_p_loc];
    B_vec = [-1/h*ones(n_p,1) 1/h*ones(n_p,1)];
    
    % Determine the number of particles in each element
    num_particles_in_e = sum(elements_particles,2);
    
    % Initialisation of element matrix and vectors
    M_e = zeros(2,2,n_e); 
    P_e = zeros(2,n_e); 
    F_grav_e = zeros(2,n_e); 
    F_int_e = zeros(2,n_e);
    
    % Determine element (lumped) matrix and vectors for active elements
    for el = 1:n_e
        if num_particles_in_e(el) > 0 
            for p = element_particle(elements_particles(el,:))
                M_e(1,1,el) = M_e(1,1,el) + mass_p(p)*N_vec(p,1);
                M_e(2,2,el) = M_e(2,2,el) + mass_p(p)*N_vec(p,2);
                P_e(:,el) = P_e(:,el) + mass_p(p)*(N_vec(p,:))'*v_p(p,s);
                F_grav_e(:,el) = F_grav_e(:,el) + (N_vec(p,:))'*f_grav_p(p);
                F_int_e(:,el) = F_int_e(:,el) + (B_vec(p,:))'*stress_p(p)*volume_p(p);
            end
            clear p
        end
    end
    clear el
    
    % Determine element vector for traction force
    F_trac_e = zeros(2,n_e); 
    el = particle_element(elements_particles(:,end));
    %f_trac_p(end) = load;
    f_trac_p(end) = tau*sin(pi*(s*t_step));
    F_trac_e(:,el) = N_vec(end,:)'*f_trac_p(end);
    clear el
    
    
    % Assemble element matrices and vectors
    M = assemble_matrix(n_n, n_e, elements_nodes_index, M_e);
    P = assemble_vector(n_n, n_e, elements_nodes_index, P_e);
    F_grav = assemble_vector(n_n, n_e, elements_nodes_index, F_grav_e);
    F_int = assemble_vector(n_n, n_e, elements_nodes_index, F_int_e);
    F_trac = assemble_vector(n_n, n_e, elements_nodes_index, F_trac_e);
    
    % Determine total force
    F = F_trac + F_grav - F_int;
    
    % Choose node to examine
    node = floor((n_e+1)/2);
    F_int_plot(1,s) = F_int(node);
    F_ext_plot(1,s) = F_grav(node);
    F_diff(1,s) = F_int_plot(1,s) - F_ext_plot(1,s);
    
    % Find inactive nodes
    non_active_nodes = find(diag(M)==0);
    
    %if langranian == 1
    %    for n = 1:n_n
    %        if ~any(non_active_nodes == n)
    %            v_n(n,s) = P(n)/M(n,n);
    %        end
    %    end
    %end
    
   
    % Calculate the nodal acceleration
    a_n = zeros(n_n,1);
    for n = 1:n_n
        if ~any(non_active_nodes == n) 
            a_n(n) = F(n)/M(n,n) ;
        end
    end
    clear n
    
    % Boundary condition
    a_n(1) = 0;
    %a_n(n_n) = 0;
    
    % Update particle velocity
    for p = 1:n_p
        v_p(p,s+1) = v_p(p,s) + t_step*N_vec(p,1)*a_n(...
            elements_nodes_index(particle_element(elements_particles...
            (:,p)),1))...
            + t_step*N_vec(p,2)*a_n(elements_nodes_index(particle_element...
            (elements_particles(:,p)),2));
    end
    clear p
    
    % Calculate updated momentum in active elements
    Mv_e = zeros(2,n_e);
    for el = 1:n_e
        if num_particles_in_e(el) > 0 
            for p = element_particle(elements_particles(el,:))
                Mv_e(:,el) = Mv_e(:,el) + mass_p(p)*(N_vec(p,:))'*v_p(p,s+1);
            end
            clear p
        end
    end
    clear el
    
    % Assemble elements vectors
    Mv = assemble_vector(n_n, n_e, elements_nodes_index, Mv_e);
    
    % Calculate the nodal velocities
    if langranian == 0
    for n = 1:n_n
        if ~any(non_active_nodes == n)
           v_n(n,s+1) = Mv(n)/M(n,n);
        end   
    end
    clear n
    end
    
    if langranian == 1
        for n = 1:n_n
            if ~any(non_active_nodes == n)
                v_n(n,s+1) = v_n(n,s) + a_n(n)*t_step;
            end
        end
    end
    
    % Boundary condition
    v_n(1,s+1) = 0;
    %v_n(n_n) = 0;

    % Calculate the nodal incremental displacement
    du_n = v_n(:,s+1)*t_step;
    
    % Boundary condition
    du_n(1) = 0;
    %du_n(n_n) = 0;

    % Update the nodal displacement
    u_n = u_n + du_n;
    
    % Determine the velocity gradient 
    %v_grad_p = zeros(n_p,1);
    for p = 1:n_p
        v_grad_p(p,1) = B_vec(p,:)*[v_n(elements_nodes_index...
            (particle_element(elements_particles(:,p)),1),s+1);
            v_n(elements_nodes_index...
            (particle_element(elements_particles(:,p)),2),s+1)];
    end
    clear p
   
    
    % Determine deformation gradient
    F_tilde_p = zeros(n_p,1);
    for p = 1:n_p
       F_tilde_p(p,1) = 1 + v_grad_p(p,1)*t_step;
    end
    clear p
   
    F_p(:,1) = ones(n_p,1);
    for p = 1:n_p
        F_p(p,s+1) = F_tilde_p(p,1).*F_p(p,s);
    end
    clear p
    
    
   % Update the mesh according to nodal velocity
   if reset == 1
       elements_nodes(:,1) = elements_nodes(:,1) + v_n(1:n_e,s+1)*t_step;
       elements_nodes(:,2) = elements_nodes(:,2) + v_n(2:n_e+1,s+1)*t_step;
   end
   
    
    % Calculate the incremental particle strain
    dstrain_p = zeros(n_p,1);
    for p = 1:n_p
         el = particle_element(elements_particles(:,p)); 
         h_temp = elements_nodes(el,2) - elements_nodes(el,1);
         B_temp = [-1/h_temp, 1/h_temp];
         dstrain_p(p) = B_temp*[du_n(elements_nodes_index...
                (particle_element(elements_particles(:,p)),1));
                du_n(elements_nodes_index...
                (particle_element(elements_particles(:,p)),2))];
    end
    clear p
    
    % Update particle stress
    stress_p = stress_p + E*dstrain_p;
    %stress_p = E.*(F_p(:,s+1)-ones(n_p,1));
    
    % Reset mesh once strain/stress have been calculated 
    if reset == 1
       elements_nodes = elements_nodes_original;  
    end 
    
    
    % Update particle position and displacement
    for p = 1:n_p
        %u_p(p,s+1) = u_p(p,s) + N_vec(p,1)*du_n(elements_nodes_index...
        %    (particle_element(elements_particles(:,p)),1))...
        %    + N_vec(p,2)*du_n(elements_nodes_index(particle_element...
        %    (elements_particles(:,p)),2));
          u_p(p,s+1) = u_p(p,s) + t_step*N_vec(p,1)*v_n(...
            elements_nodes_index(particle_element(elements_particles...
            (:,p)),1),s+1)...
            + t_step*N_vec(p,2)*v_n(elements_nodes_index(particle_element...
            (elements_particles(:,p)),2),s+1);
        
        if change_glob_pos == 1
          %  pos_p_glob(p) = pos_p_glob(p) +  N_vec(p,1)*du_n...
          %      (elements_nodes_index(particle_element...
          %      (elements_particles(:,p)),1))...
          %     + N_vec(p,2)*du_n(elements_nodes_index(particle_element...
          %      (elements_particles(:,p)),2));
            pos_p_glob(p,s+1) = pos_p_glob(p,s) + t_step*N_vec(p,1)*v_n(...
            elements_nodes_index(particle_element(elements_particles...
            (:,p)),1),s+1)...
            + t_step*N_vec(p,2)*v_n(elements_nodes_index(particle_element...
            (elements_particles(:,p)),2),s+1);
        end
    end
    clear p
   
      % Update the volume and density of the particles
   for p = 1:n_p
        volume_p(p) = (1+dstrain_p(p))*volume_p(p);
        %volume_p(p)  = det(F_p(p,s+1))*volume_p(p);
        density_p(p) = density_p(p)/(1+dstrain_p(p));
   end
   clear p 
    
    % Check if there are particles that have left their original element
    for p = 1:n_p
        el = particle_element(elements_particles(:,p));
        if pos_p_glob(p) < elements_nodes(el,1)
            fprintf('Particle crossed lower border!')
            indicator(el,s) = 0.07;
            old_elem = el;
            if old_elem > 1
                new_elem = old_elem - 1;
                elements_particles(old_elem,p) = 0;
                elements_particles(new_elem,p) = 1;              
            else
                fprintf('Particle crossed the lower boundary!')
            end
            clear old_elem new_elem
        elseif pos_p_glob(p) > elements_nodes(el,2)
            fprintf('Particle crossed upper border!')
            indicator(el+1,s) = 0.07;
            old_elem = el;
            if old_elem < n_e
                new_elem = old_elem + 1;
                elements_particles(old_elem, p) = 0;
                elements_particles(new_elem,p) = 1;
            else
                fprintf('Particle crossed the upper bundary!')
            end
            clear old_eleme new_elem
        end
    end
    clear p
    
    % Change the local positions of the particles
    for p = 1:n_p
        pos_p_loc(p) = (pos_p_glob(p) - elements_nodes(particle_element...
            (elements_particles(:,p)),1))/h;
    end
    
    clear p
    
    
    
    particles_in_element = sum(elements_particles,2);
    
    for p=1:n_e
       if particles_in_element(p) == 0
           disp('Empty cell!')
           s;
       end
    end
  
    
end