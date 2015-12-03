% Determine whether or not ULFEM_1D show logarithmic behaviour
close all;
clear all;
% Load applied on oedometer
load = [0 100 200 300 400 500 600 700 1000 10000]

% Average position of node 9, density = 1E3 Youngs modulus = 1E5
average_pos = [0.468366 0.46797 0.46758 0.46720 0.46683 0.46646 0.46610 0.46576 0.46477];



% Strain
strain1 = 0.5*ones(1,length(average_pos))- average_pos;


% Strain with large deformations
strain_l = -log(1-strain1);



figure(1)
plot(load,strain_l,'LineWidth',2)
xlabel('Load p_0 [Pa]')
ylabel('Strain \epsilon')
legend('\rho = 1E3')
title('Relation Load and Strain')

figure(2)
plot(load2,strain_l2,'LineWidth',2)
xlabel('Load p_0 [Pa]')
ylabel('Strain \epsilon')
legend('\rho = 5E3')
title('Relation Load and Strain')
