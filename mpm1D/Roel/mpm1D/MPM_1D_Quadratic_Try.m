%% MPM_1D_Quadratic_Try (uses element_particle.m, particle_element.m, assemble_matrix, assemble_vector)
%  Roel Tielen
%  November 19th 2016
%  TU Delft
%  Based on the initial code of Lisa Wobbes. 

function[u_n, v_p, v_n, mass_p,u_p,E_kin,E_pot,E_grav,E_trac, F_diff, F_int_plot,F_ext_plot, stress_particle,indicator] = MPM_1D_Quadratic_Try(density, E, g, load, H, n_e, h, n_ep,...
    pos_p_glob,  pos_p_loc, t_step, n_time_steps, total_time, u_p, v_p, stress_p, u_n,...
    v_n, both_ends_fixed, change_glob_pos,change_loc_pos, E_kin,E_pot,E_grav,E_trac,lagranian, volume_update,lumped, ULFEM,momentum, alpha) 
 % INPUT: density, Young's modulus, gravitational acceleration, load,
 % height/length, number of elements, element size, number of particles per
 % element, global particle position, local particle positions, time step, 
 % number of time steps, total time, initial particle displacement, initial 
 % particle velocity, initial particle stress, initial nodal displacement, 
 % initial nodal velocity, bounday conditions, update global positions, update local positions, 
 % kinetic energy, potential energy, gravitational energy, traction energy, MPM algorithm, volume updating)
 % OUTPUT: nodal displacement, particle velocity, nodal velocity, particle
 % mass, particle displacement, kinetic energy, potential energy,
 % gravitational energy, traction energy, force difference, internal
 % force,external force, particle stress, grid crossing indicator.



%% Mesh and element connections
% Define nodes and mesh
n_n = 2*n_e+1; 
mesh  = 0:h:H; 

% Define elements and index of elements
elements_nodes = zeros(n_e,3);
elements_nodes_index = zeros(n_e,3);
for el = 1:n_e
    elements_nodes(el,:) = [mesh(el) (mesh(el)+mesh(el+1))/2 mesh(el+1)]; 
    elements_nodes_index(el,:) = [2*el-1 2*el 2*el+1]; 
end
clear el


%% Particle initialization
% Number of particles
n_p = n_e*n_ep; 

% Initial particle location
elements_particles = zeros(n_e,n_p);

 for i=1:n_e
    for j = (i-1)*n_ep+1:i*n_ep
      elements_particles(i,j) = 1;  
    end
 end

% Initial particle properties
volume_p = (1/n_ep)*h*ones(n_p,1); 
density_p = density*ones(n_p,1); 
mass_p = volume_p.*density_p; 
f_grav_p = mass_p*g; 
f_trac_p = zeros(n_p,1); 
f_trac_p(end) = load; 
v_n = zeros(n_n,n_time_steps);

% Vectors to store information
F_int_plot = zeros(1,n_time_steps);
F_ext_plot = zeros(1,n_time_steps);
F_diff = zeros(1,n_time_steps);
stress_particle = zeros(1,n_time_steps);
indicator = zeros(n_n,n_time_steps);

%% Time integration
for s = 1:n_time_steps-1
    
    % Basis functions and derivatives in material points
    N_vec = [2*(ones(n_p,1)/2 - pos_p_loc).*(ones(n_p,1) - pos_p_loc) 4.*pos_p_loc.*(ones(n_p,1)-pos_p_loc)    -2*pos_p_loc.*(ones(n_p,1)/2 - pos_p_loc)]; 
    B_vec = [4*pos_p_loc- 3*ones(n_p,1)  4*ones(n_p,1) - 8*pos_p_loc 4*pos_p_loc - ones(n_p,1) ];
    
    % Define load on last particle
    f_trac_p(end) = load; 
    
    % Determine the number of particles in each element
    num_particles_in_e = sum(elements_particles,2);
  
    % Initialisation of element matrix and vectors
    M_e = zeros(3,3,n_e); 
    P_e = zeros(3,n_e); 
    F_grav_e = zeros(3,n_e); 
    F_int_e = zeros(3,n_e); 
    
    % Determine element (lumped) matrix and vectors for active elements
    for el = 1:n_e
        if num_particles_in_e(el) > 0 
            for p = element_particle(elements_particles(el,:))
                if ULFEM == 1
                    h = elements_nodes(el,3) - elements_nodes(el,1);
                    if lumped == 1                       
                        M_e(1,1,el) = M_e(1,1,el) + h*density_p(p)*N_vec(p,1);
                        M_e(2,2,el) = M_e(2,2,el) + h*density_p(p)*N_vec(p,2);
                        M_e(3,3,el) = M_e(3,3,el) + h*density_p(p)*N_vec(p,3);
                    else
                        M_e(1,1,el) = M_e(1,1,el) + h*density_p(p)*N_vec(p,1)*N_vec(p,1);
                        M_e(1,2,el) = M_e(1,2,el) + h*density_p(p)*N_vec(p,1)*N_vec(p,2);
                        M_e(1,3,el) = M_e(1,3,el) + h*density_p(p)*N_vec(p,1)*N_vec(p,3);
                        M_e(2,1,el) = M_e(2,1,el) + h*density_p(p)*N_vec(p,2)*N_vec(p,1);
                        M_e(2,2,el) = M_e(2,2,el) + h*density_p(p)*N_vec(p,2)*N_vec(p,2);
                        M_e(2,3,el) = M_e(2,3,el) + h*density_p(p)*N_vec(p,2)*N_vec(p,3);
                        M_e(3,1,el) = M_e(3,1,el) + h*density_p(p)*N_vec(p,3)*N_vec(p,1);
                        M_e(3,2,el) = M_e(3,2,el) + h*density_p(p)*N_vec(p,3)*N_vec(p,2);
                        M_e(3,3,el) = M_e(3,3,el) + h*density_p(p)*N_vec(p,3)*N_vec(p,3);
                    end
                    P_e(:,el)       = P_e(:,el) + h*density_p(p)*(N_vec(p,:))'*v_p(p,s); 
                    F_grav_e(:,el)  = F_grav_e(:,el) + (N_vec(p,:))'*h*density_p(p)*g;
                    F_int_e(:,el)   = F_int_e(:,el) + (B_vec(p,:))'*stress_p(p)*h;    
                else    
                    if lumped == 1
                        M_e(1,1,el) = M_e(1,1,el) + mass_p(p)*N_vec(p,1);
                        M_e(2,2,el) = M_e(2,2,el) + mass_p(p)*N_vec(p,2);
                        M_e(3,3,el) = M_e(3,3,el) + mass_p(p)*N_vec(p,3);
                    else
                        M_e(1,1,el) = M_e(1,1,el) + mass_p(p)*N_vec(p,1)*N_vec(p,1);
                        M_e(1,2,el) = M_e(1,2,el) + mass_p(p)*N_vec(p,1)*N_vec(p,2);
                        M_e(1,3,el) = M_e(1,3,el) + mass_p(p)*N_vec(p,1)*N_vec(p,3);
                        M_e(2,1,el) = M_e(2,1,el) + mass_p(p)*N_vec(p,2)*N_vec(p,1);
                        M_e(2,2,el) = M_e(2,2,el) + mass_p(p)*N_vec(p,2)*N_vec(p,2);
                        M_e(2,3,el) = M_e(2,3,el) + mass_p(p)*N_vec(p,2)*N_vec(p,3);
                        M_e(3,1,el) = M_e(3,1,el) + mass_p(p)*N_vec(p,3)*N_vec(p,1);
                        M_e(3,2,el) = M_e(3,2,el) + mass_p(p)*N_vec(p,3)*N_vec(p,2);
                        M_e(3,3,el) = M_e(3,3,el) + mass_p(p)*N_vec(p,3)*N_vec(p,3);                        
                    end
                    P_e(:,el)       = P_e(:,el) + mass_p(p)*(N_vec(p,:))'*v_p(p,s); 
                    F_grav_e(:,el)  = F_grav_e(:,el) + (N_vec(p,:))'*f_grav_p(p);
                    F_int_e(:,el)   = F_int_e(:,el) + (B_vec(p,:))'*stress_p(p)*volume_p(p);                
                end
            end
            clear p
        end
    end
    clear el
  
    % Determine element vector for traction force
    F_trac_e = zeros(3,n_e); 
    el = particle_element(elements_particles(:,end)); 
    F_trac_e(:,el) = (N_vec(end,:))'*f_trac_p(end);
    clear el
    
    % Assemble element matrices and vectors
    M = assemble_matrix_quadratic(n_n, n_e, elements_nodes_index, M_e);
    F_grav = assemble_vector_quadratic(n_n, n_e, elements_nodes_index, F_grav_e);
    F_int = assemble_vector_quadratic(n_n, n_e, elements_nodes_index, F_int_e);
    P = assemble_vector_quadratic(n_n, n_e, elements_nodes_index, P_e);
    F_trac = assemble_vector_quadratic(n_n, n_e, elements_nodes_index, F_trac_e);
    
    % Determine total force
    F = F_trac + F_grav - F_int;
    
    % Choose node to examine
    node = floor(n_n/2);
    F_int_plot(1,s) = F_int(node);
    F_ext_plot(1,s) = F_grav(node);
    F_diff(1,s) = F_int_plot(1,s) - F_ext_plot(1,s);
    
    % Apply boundary condition on momentum
    if both_ends_fixed == 1 
        P(1) = 0;
        P(n_n) = 0;
    else 
        P(1) = 0;
    end
    
    % Find inactive nodes
    non_active_nodes = find(diag(M)==0);
   
    % Determine nodal velocity
    if lagranian == 1
        if momentum == 1
        for n = 1:n_n
            if ~any(non_active_nodes == n)
                v_n(n,s) = P(n)/M(n,n);
            end
        end
        end
    end
    clear n 
    
    % Calculate the nodal acceleration
    a_n = zeros(n_n,1);
    for n = 1:n_n
        if lumped == 1 
            if ~any(non_active_nodes == n) 
                a_n(n) = F(n)/M(n,n) ;
            end
        else
            a_n = M\F;
        end
    end
    clear n
    
    % Apply boundsry condition on nodal acceleration
    if both_ends_fixed == 1
        a_n(1,1) = 0;
        a_n(n_n) = 0;
    else 
        a_n(1,1) = 0;
    end
    
    % Update particle velocity
    for p = 1:n_p
        v_p(p,s+1) = v_p(p,s) + t_step*N_vec(p,1)*a_n(...
            elements_nodes_index(particle_element(elements_particles...
            (:,p)),1))+ t_step*N_vec(p,2)*a_n(elements_nodes_index...
            (particle_element(elements_particles(:,p)),2)) + t_step*N_vec(p,3)...
            *a_n(elements_nodes_index(particle_element(elements_particles...
            (:,p)),3));
    end
    clear p
  
    % Calculate updated momentum in active elements
    Mv_e = zeros(3,n_e);
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
    Mv = assemble_vector_quadratic(n_n, n_e, elements_nodes_index, Mv_e);
    
    % Apply boundary condition on updated momentum
    if both_ends_fixed == 1
        Mv(1)=0;
        Mv(n) = 0;
    else
        Mv(1) = 0;
    end
    
    % Calculate the nodal velocities at new time step
    if lagranian == 0
    for n = 1:n_n
        if ~any(non_active_nodes == n)
           v_n(n,s+1) = Mv(n)/M(n,n);
        end   
    end
    clear n
    end
    
    if lagranian == 1
        for n = 1:n_n
                v_n(n,s+1) = v_n(n,s) + a_n(n,1)*t_step;
        end
    end
    clear n 
    
    % Apply boundary condition on nodal velocity
    if both_ends_fixed == 1
        v_n(1,:) = 0;
        v_n(n_n) = 0;
    else
        v_n(1,:) = 0;
    end
    

    % Calculate the nodal incremental displacement
    du_n = v_n(:,s+1)*t_step;

    % Update the nodal displacement
    u_n(:,s+1) = u_n(:,s) + du_n;
    
    % Apply boundary condition on nodal displacement
    if both_ends_fixed == 1
        u_n(1,s+1) = 0;
        u_n(n_n,s+1) = 0;
    else
        u_n(1,s+1) = 0;
    end
    
    
    % Calculate the incremental particle strain
    dstrain_p = zeros(n_p,1);
    for p = 1:n_p
        dstrain_p(p) = B_vec(p,:)*[du_n(elements_nodes_index...
            (particle_element(elements_particles(:,p)),1));
            du_n(elements_nodes_index...
            (particle_element(elements_particles(:,p)),2));...
            du_n(elements_nodes_index...
            (particle_element(elements_particles(:,p)),3))];
    end
    clear p
    
    % Update particle stress
    stress_p = stress_p + E*dstrain_p;
    
    stress_particle(1,s) = stress_p(1);
   
    % Update particle position and displacement
    for p = 1:n_p
            u_p(p,s+1) = u_p(p,s) + N_vec(p,1)*du_n(elements_nodes_index...
                (particle_element(elements_particles(:,p)),1))...
                + N_vec(p,2)*du_n(elements_nodes_index(particle_element...
                (elements_particles(:,p)),2)) + N_vec(p,3)*du_n(elements_nodes_index...
                (particle_element(elements_particles(:,p)),3));
        if change_glob_pos == 1
            pos_p_glob(p) = pos_p_glob(p) +  N_vec(p,1)*du_n...
                (elements_nodes_index(particle_element...
                (elements_particles(:,p)),1))...
                + N_vec(p,2)*du_n(elements_nodes_index(particle_element...
                (elements_particles(:,p)),2)) + N_vec(p,3)*du_n...
                (elements_nodes_index(particle_element...
                (elements_particles(:,p)),3));
        end
    end
    clear p
    
   % Update the volume and density of the particles
   if volume_update == 1
    for p = 1:n_p
        volume_p(p) = (1+dstrain_p(p))*volume_p(p);
        density_p(p) = density_p(p)/(1+dstrain_p(p));
    end
    clear p 
   end
   
   % Update nodal positions according to nodal displacement
   if ULFEM == 1
       elements_nodes(:,1) = elements_nodes(:,1) + v_n(1:n_e,s+1)*t_step;
       elements_nodes(:,2) = elements_nodes(:,2) + v_n(2:n_e+1,s+1)*t_step;
   end
   
   % Check if there are particles that have left their original element
   for p = 1:n_p
       el = particle_element(elements_particles(:,p));
       if pos_p_glob(p) < elements_nodes(el,1)
           fprintf('Particle crossed lower border!')
           indicator(el,s) = 350;
           old_elem = el;
           if old_elem > 1
                new_elem = old_elem - 1;
                elements_particles(old_elem,p) = 0;
                elements_particles(new_elem,p) = 1;
                time_lower_crossing = s*t_step
            else
                fprintf('Particle crossed the lower boundary!')
            end
            clear old_elem new_elem
        elseif pos_p_glob(p) > elements_nodes(el,2)
            fprintf('Particle crossed upper border!')
            indicator(el,s) = 350;
            old_elem = el;
            if old_elem < n_e
                new_elem = old_elem + 1;
                elements_particles(old_elem, p) = 0;
                elements_particles(new_elem,p) = 1;
                time_upper_crossing = s*t_step
            else
                fprintf('Particle crossed the upper bundary!')
            end
            clear old_eleme new_elem
        end
    end
    clear p
   
    % Change the local positions of the particles
    if change_loc_pos == 1
        for p = 1:n_p
            pos_p_loc(p) = (pos_p_glob(p) - elements_nodes(particle_element...
            (elements_particles(:,p)),1))/h;
        end
    end
    clear p
    
    % Check whether an empty cells occurs and print time
    particles_in_element = sum(elements_particles,2);
    
    for el=1:n_e
       if particles_in_element(el) == 0
           disp('Empty cell!');
           time_empty_cell = s*t_step
           number_elements_in_particle = sum(elements_particles,2);
       end
    end
    clear el
    
    % Determine the kinetic and potential energy of the system
    for n = 1:n_n
       E_kin(1,s+1) = E_kin(1,s+1) + 0.5*M(n,n)*v_n(n,s+1)*v_n(n,s+1); 
       E_pot(1,s+1) = E_pot(1,s+1) + -0.5*F_int(n)*u_n(n,s+1);
       E_trac(1,s+1) = E_trac(1,s+1) + F_trac(n)*u_n(n,s+1);
       E_grav(1,s+1) = E_grav(1,s+1) + F_grav(n)*u_n(n,s+1);
    end
end

