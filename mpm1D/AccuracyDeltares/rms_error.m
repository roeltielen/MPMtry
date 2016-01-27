function [] = rms_error()
%% Provide required info
% initial particle positions
x0 = read_position('x_init_fem_128.txt');
x0 = x0*10^(-8);

%x01 = 
% other parameters
L = 25;
rho = 1;
E = 100;
v0 = 0.1;
T = 0.5;
t_step = 1E-3;
n = T/t_step;

%% Compute the exact solution
x_exact = exact_solution(x0, n, t_step, v0, E, rho, L);

%% Read the numerical solution
x_num = read_position('x_end_fem_128.txt');
x_num = x_num*10^(-8);
%% Compute the RMS error
number_particles = size(x0,1);
Error = norm((x_exact-x_num))/sqrt(number_particles)

