% Determine the analytical solution of the strain as function of the load
close all;
clear all;

%% Parameters
% Define parameters
E = 1E5;
p_0 = 0:1:30000;


%% Analytical solution
% Calculate strain for FEM and ULFEM
strainFEM_exact = p_0/E;
strainULFEM_exact = log(ones(1,length(p_0)) + p_0/E);

%% Numerical solution
% Average position of node 9, density = 1E3, Youngs modulus = 1E5 no gravity!
load = [0 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 11000 12000 13000 14000 15000];
average_FEM = [0.5 0.4950 0.4900 0.4850 0.4800 0.4750 0.4700 0.4650 0.4600 0.4550 0.4500 0.4450 0.4400 0.4350 0.4300 0.4250 ]; 

% Determine strain with FEM code and ULFEM code
strainFEM = (0.5*ones(1,length(average_FEM))- average_FEM)/0.5;
strainULFEM = [0 0.0099 0.0198 0.0296 0.0392 0.0488 0.0583 0.0677 0.0770 0.0862 0.0954 0.1044 0.1134 0.1223 0.1312 0.1399];

% Determine strain with MPM-ULFEM, t = 2.5 s, delta t = 1E-4 s, 16 elements
load2 =           [0  5000     10000   15000   20000     25000   30000];
strainULFEM_MPM = [0  0.04798  0.0943  0.1382  0.181625  0.2226  0.2626 ];

% Strain MPM with small deformation!
strainMPM =       [0  ];

% Determine strain with MPM-FEM, t = 2.5 s, delta t = 1E-3
average_FEM_MPM = [0.5 0.4750 0.4500 0.4250 0.4000 0.3750 0.3500 ]; 
strainFEM_MPM = (0.5*ones(1,length(average_FEM_MPM))- average_FEM_MPM)/0.5;

% Determine strain with ULFEM code
strain3 = -log(1-strainFEM);




%% Plot 
% Plot solution
plot(strainFEM_exact,p_0,'--r','LineWidth',1.25)
hold on
plot(strainULFEM_exact,p_0,'--b','LineWidth',1.25)
%plot(strainFEM,load,'*')
%plot(strainULFEM,load,'*r')
plot(strainFEM_MPM,load2,'*')
plot(strainULFEM_MPM,load2,'*r')
ylabel('Load p_0 [Pa]','Fontsize',12)
xlabel('Strain \epsilon','Fontsize',12)
legend('Exact FEM','Exact ULFEM','FEM','ULFEM','Location','northwest')
title('Relation Load and Strain','Fontsize',12)
