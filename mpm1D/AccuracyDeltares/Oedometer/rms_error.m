function [] = rms_error()
%% Provide required info
% initial particle positions
x0 =  read_position('x_init_fem_64.txt');
x0 = x0*10^(-8);

%x01 = 
% other parameters
H = 25;
rho = 1000;
E = 1E5;
g = 9.81;
p0 = 0;
T = 0.5;
%t_step = 1E-3;
%n = T/t_step;

%% Compute the exact solution
number_nodes = size(x0,1);
x_exact = zeros(number_nodes,1);
for i = 1:number_nodes 
    x_exact(i) = exact_solution(rho, E, p0, g, H, x0(i),T);
end

%% Read the numerical solution
x_num = read_position('x_end_fem_64gl.txt');
x_num = x_num*10^(-8);
%% Compute the RMS error
%number_particles = size(x0,1);
Error = norm((x_exact-x_num))/sqrt(number_nodes)
Rel_error = Error/(norm(x_num)/sqrt(number_nodes))

