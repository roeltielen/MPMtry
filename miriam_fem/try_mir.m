function [num_sol, exact_sol, M, time_step, dt] = try_mir()

clear all
clc

%% Get input
n_e = 4;
h = 1/n_e;
H = 1;


nod  = 0:h:H; %global mesh
nod = nod';
n_nod = length(nod);
elements = zeros(n_e,2);
elm = zeros(n_e,2);
for i = 1:n_e
   elements(i,:) = [nod(i) nod(i+1)]; 
   elm(i,:) = [i i+1]; %contains the indices of the
    %nodes comprising the elements
end
n_elm = length(elm);

rho = 1E3;    % Density Solid
E = 1E5;      % Young's Modulus
nu = 0;     % Poisson Ratio

u = zeros(n_e+1,1);
v = zeros(n_e+1,1);
s = zeros(n_e,1);

fix = [1 1];
n_fix = 1;

load = [n_e+1 -5E3];
n_load = 1;
p = -5E3;

dt = 1E-4;   % Size of time step
m = 600;    % Number of time steps
m = 100*m;
a = 0;    % Damping factor

%% Define other variables

g = 0;        % Gravitational acceleration

%% Define shape function for elements
% First order interpolation functions
% Single Gauss Point in normalized linear element

N = 0.5;
LN = [-1,1];

%% Define shape function for boundary elements
% First order interpolation functions
% Single Gauss Point in normalized point element

N_b = 1;

%% Define lumped-mass matrix (time-independent)

M = zeros(n_nod,1);

% for-loop over elements
for i=1:n_elm
    getnod = elm(i,:);
    J = nod(getnod(2)) - nod(getnod(1));
    M(getnod(1)) = M(getnod(1)) + N * rho * J;
    M(getnod(2)) = M(getnod(2)) + N * rho * J;
end

%% Define traction force vector (time-independent)

F_trac = zeros(n_nod,1);

% for-loop over boundary elements
for i=1:n_load
    getload = load(i,:);
    F_trac(getload(1)) = F_trac(getload(1)) + N_b * getload(2);
end

%% Define gravity force vector (time-independent)

F_grav = - M .* g;


%% TIME INTEGRATION

saveu = zeros(n_nod,m+1);
saveu(:,1) = nod;
savev = zeros(n_nod,m+1);
saves = zeros(n_elm,m+1);

for t = 1:m
    %% Define internal force vector (time-dependent)

    F_int = zeros(n_nod,1);

    % for-loop over elements
    for i=1:n_elm
        getnod = elm(i,:);
        F_int(getnod(1)) = F_int(getnod(1)) + LN(1)*s(i);
        F_int(getnod(2)) = F_int(getnod(2)) + LN(2)*s(i);
    end
    
    %% Update velocity, stress and displacement
    
    % update velocity
    v = v + M .\ (F_trac - F_int + F_grav) .* dt;
    
    % fixities
    for i = 1:n_fix
        if(fix(i,2)==1)
            v(fix(i,1)) = 0;
        end
    end
    
    % calculate strain rate
    de = zeros(n_elm,1);
    for i=1:n_elm
        getnod = elm(i,:);
        J = nod(getnod(2)) - nod(getnod(1));
        de(i) = (LN(1) * v(getnod(1)) + LN(2) * v(getnod(2)))/J;
    end
    
    % update stress
    s = s + dt * E * de;
    
    % update displacement
    u = u + dt * v;    
    
    %saveu(:,t+1) = nod + u;
    saveu(:,t+1) = u;
    savev(:,t+1) = v;
    saves(:,t+1) = s;
end
format long
 saveu(:,1:20)

time_step =  510;
%t =  0 : dt : m*dt;
T_ime = time_step*dt;
time_step =  time_step+1;
for node = 1:n_e + 1
    [position_exactT(node,:),displacement_exactT(node,:),...
        velocity_exactT(node,:)] = exact_solution...
        (rho,E,p,g,H,nod(node), T_ime);
end
clear nodexact_sol
exact_sol = displacement_exactT;
num_sol = saveu(:, time_step);

errnrm = compute_error_norm(num_sol,exact_sol, diag(M));
fprintf('For time t = %e\n', T_ime)
fprintf('The error norm = %e\n', errnrm)



