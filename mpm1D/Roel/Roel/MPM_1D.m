function[u_p, v_p] = MPM_1D(density, E, g, load, H, n_e, h, n_ep,...
    pos_p_glob, t_step, n_time_steps, total_time, u_p, v_p, stress_p, u_n,...
    v_n, both_ends_fixed, change_glob_pos, pos_p_loc) 

%% Mesh and element connections
n_n = n_e+1; % number of nodes
%..Generate global mesh, element connections
mesh  = 0:h:H; %global mesh
mesh_updated = mesh'; %required to store new global positions of the nodes
elements_nodes = zeros(n_e,2);
elements_nodes_index = zeros(n_e,2);
for el = 1:n_e
    elements_nodes(el,:) = [mesh(el) mesh(el+1)]; %contains the coordinates
    %of the nodes comprising the element
    elements_nodes_index(el,:) = [el el+1]; %contains the indices of the
    %nodes comprising the elements
end
clear el

%% Particle initialization
n_p = n_e*n_ep; %total number of particles

%%..Initial particle locations
elements_particles(:,:,1) = eye(n_e,n_p); %matrix containing the
%information on the location of the particles within elements;
%initially every element contains exactly one particle
%pos_p_loc = ones(n_p,1)*0.5; % initial local positions of particles
%pos_p_glob = (mesh(1:end-1))' + h*pos_p_loc; % initial global position

%%..Initial particle properties
volume_p = (1/n_ep)*h*ones(n_p,1); %particle volume
density_p = density*ones(n_p,1); %particle density
mass_p = volume_p.*density_p; %particle mass
f_grav_p = mass_p*g; %gravitational force acting on particle
f_trac_p = zeros(n_p,1); %traction force at particles
f_trac_p(end) = load; %only the top particle is traction containing


%% Time integration
for s = 1:n_time_steps-1
    
    f_trac_p(end) = load*s*t_step/total_time; %only the top particle is traction containing
    
    %..Interpolation function and strain-displacement vectors (elementwise)
    N_vec = [ones(n_p,1) - pos_p_loc pos_p_loc];
    B_vec = [-1/h*ones(n_p,1) 1/h*ones(n_p,1)];
    
    %..Compute the number of particles in element i
    num_particles_in_e = sum(elements_particles,2);
    
    %..Map particle information to nodes
    M_e = zeros(2,2,n_e); %initialize element mass matrices
    %     P_e = zeros(2,n_e); %initialize element momentum vectors
    F_grav_e = zeros(2,n_e); %initialize element gravity force vectors
    F_int_e = zeros(2,n_e); %initialize internal force vectors
    for el = 1:n_e
        if num_particles_in_e(el) > 0 %select only active elements
            for p = element_particle(elements_particles(el,:))
                M_e(1,1,el) = M_e(1,1,el) + mass_p(p)*N_vec(p,1);
                M_e(2,2,el) = M_e(2,2,el) + mass_p(p)*N_vec(p,2);
                
                %                 P_e(:,el) = P_e(:,el) + mass_p(p)*(N_vec(p,:))'*v_p(p);
                
                F_grav_e(:,el) = F_grav_e(:,el) + (N_vec(p,:))'...
                    *f_grav_p(p);
                
                F_int_e(:,el) = F_int_e(:,el) + (B_vec(p,:))'...
                    *stress_p(p)*volume_p(p);
            end
            clear p
            
        end
    end
    clear el
    
    
        F_trac_e = zeros(2,n_e); %initialize traction force vectors
        el = particle_element(elements_particles(:,end)); %element containing
        %traction carrying particle
        F_trac_e(:,el) = (N_vec(end,:))'*f_trac_p(end);
        clear el
    
    %..Assemblage procedure
    M = assemble_matrix(n_n, n_e, elements_nodes_index, M_e);
    %     P = assemble_vector(n_n, n_e, elements_nodes_index, P_e);
    F_grav = assemble_vector(n_n, n_e, elements_nodes_index, F_grav_e);
    F_int = assemble_vector(n_n, n_e, elements_nodes_index, F_int_e);
    %     F_trac = assemble_vector(n_n, n_e, elements_nodes_index, F_trac_e);
    
    %..Compute the total force
    %     F = F_trac + F_grav - F_int;
    F = F_grav - F_int;
    
    %..Find inactive nodes
    non_active_nodes = find(diag(M)==0);
    
    %..Calculate the nodal acceleration
    a_n = zeros(n_n,1);
    for n = 1:n_n
        if ~any(non_active_nodes == n) %M(n,n) > 0.1*sum(diag(M))
            a_n(n) = F(n)/M(n,n) ;
%         else
%             a_n(n) = 0.5;
        end
    end
    clear n
    
    %..Boundary condition
    a_n(1) = 0;
    a_n(n_n) = 0;
    
    %..Update particle velocity
    for p = 1:n_p
        v_p(p,s+1) = v_p(p,s) + t_step*N_vec(p,1)*a_n(...
            elements_nodes_index(particle_element(elements_particles...
            (:,p)),1))...
            + t_step*N_vec(p,2)*a_n(elements_nodes_index(particle_element...
            (elements_particles(:,p)),2));
    end
    clear p
    
    %..Calculate updated momentum
    Mv_e = zeros(2,n_e);
    for el = 1:n_e
        if num_particles_in_e(el) > 0 %select only active elements
            for p = element_particle(elements_particles(el,:))
                
                Mv_e(:,el) = Mv_e(:,el) + mass_p(p)*(N_vec(p,:))'*v_p(p,s+1);
                
            end
            clear p
            
        end
    end
    clear el
    %..Assemble
    Mv = assemble_vector(n_n, n_e, elements_nodes_index, Mv_e);
    
    %..Obtain the nodal velocities
    v_n = zeros(n_n,1);
    for n = 1:n_n
        if ~any(non_active_nodes == n)
            v_n(n) = Mv(n)/M(n,n);
%         else
%             v_n(n) = a_n(n)*t_step;
        end
    end
    clear n
    v_n(1) = 0;
    v_n(n_n) = 0;
    %     %..Boundary condition
    %     v_n(1) = 0;
    %
    
    %..Calculate the nodal incremental displacement
    du_n = v_n*t_step;
    
    %..Boundary condition
    du_n(1) = 0;
    du_n(n_n) = 0;
    %     for n = non_active_nodes
    %         du_n(n) = 0;
    %     end
    %     clear n
    
    %..Update the nodal displacement
    u_n = u_n + du_n;
    
    %     for n = non_active_nodes
    %         u_n(n) = 0;
    %     end
    %     clear n
    
    %%..New: Update vector B according to new nodal positions
    mesh_updated = mesh_updated + du_n; %update the mesh
    % Compute the new distance between the nodes
    h_updated = zeros(n_p,1);
    for p = 1:n_p
        h_updated(p) = mesh_updated(elements_nodes_index(particle_element...
            (elements_particles(:,p)),2))-...
        mesh_updated(elements_nodes_index(particle_element...
            (elements_particles(:,p)),1));
    end
    clear p
    % Update B_vec
    B_vec_new(:,1) = -1./h_updated;
    B_vec_new(:,2) = 1/h_updated;    
    
    %..Calculate the incremental particle strain
    dstrain_p = zeros(n_p,1);
    for p = 1:n_p
        dstrain_p(p) = B_vec(p,:)*[du_n(elements_nodes_index...
            (particle_element(elements_particles(:,p)),1));
            du_n(elements_nodes_index...
            (particle_element(elements_particles(:,p)),2))];
    end
    clear p
    
    %..Update particle stress
    stress_p = stress_p + E*dstrain_p;
    
    %..Update particle positiion and displacement
    for p = 1:n_p
        u_p(p,s+1) = u_p(p,s) + N_vec(p,1)*du_n(elements_nodes_index...
            (particle_element(elements_particles(:,p)),1))...
            + N_vec(p,2)*du_n(elements_nodes_index(particle_element...
            (elements_particles(:,p)),2));
        if change_glob_pos == 1
            pos_p_glob(p) = pos_p_glob(p) +  N_vec(p,1)*du_n...
                (elements_nodes_index(particle_element...
                (elements_particles(:,p)),1))...
                + N_vec(p,2)*du_n(elements_nodes_index(particle_element...
                (elements_particles(:,p)),2));
        end
    end
    clear p
    
    %%..Book-keeping
    %%..Check if there are particles that have left their original element
    for p = 1:n_p
        el = particle_element(elements_particles(:,p));
        if pos_p_glob(p) < elements_nodes(el,1)
            fprintf('Particle crossed lower border!')
            %             s
            %             p
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
            %             s
            %             p
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
    
    %%..Change the local positions of the particles
    for p = 1:n_p
        pos_p_loc(p) = (pos_p_glob(p) - elements_nodes(particle_element...
            (elements_particles(:,p)),1))/h;
    end
    clear p
    
end