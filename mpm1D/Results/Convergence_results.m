%% File in which the error of the MPM 1D solution of the vibrating string solution is stored.
% Input
% density = 100;
% Youngs_modulus = 100;
% gravitational_acceleration = 0; 
% load = 0.3;
% height = 1; 

close all
clear all

% 3 particle per elements Error at time t = 2, delta t = 0.01. Length of simulation 8 sec.
el = [10 20 40 60 80 100 120 140 160 ];
er = [5.466e-4 4.899e-4 4.394e-04 4.182e-04 4.050e-04 3.960e-04 3.727e-04 3.440e-04 3.562e-4];

% 3 particles per element Error at time t = 1, delta t = 0.001. Length of simulation 4 sec.
el2 = [10 20 40 80 160 ]
er2 = [6.9285e-04  5.6973e-04 5.2487e-04  4.9765e-04  4.8114-04 ]; 

% Error same as er2, but when bar is equally discretized with 3.5 PPC. 
el2_CFL = [10 20 40 80 160 320] 
er2_uniform = [6.26e-04 5.3900e-04 5.0548e-04 4.8445e-04 4.749e-04];

% With CFL 0.8
er_CFL = [2.8990e-04 2.8910e-04 3.7824e-04 4.1548e-04 4.4091e-04 4.707e-04];

% Explicit forward-Euler, CFL 0.5,3.5 PPC Error at time t = 1 s. Length of simulation 4 sec. 
er_CFL2 = [ 3.1336e-04 3.9332e-04 4.5365e-04 4.8925e-04 5.0210e-04  5.0850e-04];

% Explicit forward-Euler, CFL 0.5, +- 3 PPC Error at time t = 1 s. Length of simulation 4 sec. 
er_CFL3 = [ 3.5290e-04 3.8263e-04  4.4487e-04 4.7482e-04 4.8613e-04  4.9266e-04];

figure(1)
loglog(el,er,'-*','LineWidth',2)
xlabel('Number of elements','Fontsize',12)
ylabel('RMS Error','Fontsize',12)
axis([10 10e2 10e-8 10e-3])
title('Spatial convergence MPM','Fontsize',12)
legend('MPM')

figure(2)
loglog(el2,er2_uniform,'-*','LineWidth',2)
xlabel('Number of elements','Fontsize',12)
ylabel('RMS Error','Fontsize',12)
axis([10 10e2 10e-8 10e-3])
title('Spatial convergence MPM','Fontsize',12)
legend('MPM')

figure(3)
loglog(el2_CFL,er_CFL,'-*','LineWidth',2)
xlabel('Number of elements','Fontsize',12)
ylabel('RMS Error','Fontsize',12)
axis([10 10e2 10e-8 10e-3])
title('Spatial convergence MPM','Fontsize',12)
legend('MPM')

figure(4)
loglog(el2_CFL,er_CFL3,'-*','LineWidth',2)
xlabel('Number of elements','Fontsize',12)
ylabel('RMS Error','Fontsize',12)
axis([10 10e2 10e-8 10e-3])
title('Spatial convergence MPM','Fontsize',12)
legend('MPM')

%% File in which the error of the MPM 1D solution of the vibrating string solution is stored.
%% Input 
% density = 1;
% Youngs_modulus = 100;
% gravitational_acceleration = 0; 
% load = 0;
% height = 25;
% Number of elements
el = [20 40 60 80 100 120 140 160 180 200 220 240 280 ]

% 4 particles per element Error at time t = 0.5, delta t = 0.01. Length of simulation 5 sec
er4 =      [1.4546e-04  3.5744e-05  1.5405e-05  8.2875e-06  4.9977e-06   5.4264e-06 1.0289e-05 1.6661e-05 2.3662e-05 2.8483e-05 3.7626e-05 4.7159e-05 6.9502e-05 ];
vol_er4 =  [1.4555e-04  3.6251e-05  1.6606e-05  1.0384E-05  8.0279E-06   8.2979E-06 1.1893E-05 1.7469e-05 2.4064e-05 2.8713e-05 3.7743e-05 4.7219e-05 6.9340e-05 ];       

% 6 particles per element Error at time t = 0.5, delta t = 0.01. Length of simulation 5 sec
vol_er6 =  [1.3905e-04  3.4655e-05  1.5947e-05  1.0079e-05  1.1209e-05   1.5758e-05 2.2981e-05 2.8203e-05 3.3774e-05 3.7937e-05 4.2115e-05 4.6099e-05 5.4192e-05];

% 8 particles per element Error at time t = 0.5, delta t = 0.01. Length of simulation 5 sec
vol_er8 =  [1.3550e-04  3.3786e-05  1.5590e-05  1.2908e-05  1.7222e-05   2.3460e-05 2.8197e-05 3.2039e-05 3.5176e-05 3.5813e-05 3.7102e-05 3.6440e-05 4.0484e-05];   

% Particle boundary crossing starts at 80 elements 

figure(5)
semilogy(el,er4,'-*','LineWidth',2)
hold on
semilogy(el,vol_er4,'-*r','LineWidth',2)
xlabel('Number of elements','Fontsize',12)
ylabel('RMS Error','Fontsize',12)
legend('4 PPC','4 PPC Volume')
title('Effect volume updating on RMS error','Fontsize',12)

figure(6)
semilogy(el,er4,'-*r','LineWidth',2)
hold on
semilogy(el,vol_er6,'-*','LineWidth',2)
semilogy(el,vol_er8,'-*m','LineWidth',2)
xlabel('Number of elements','Fontsize',12)
ylabel('RMS Error','Fontsize',12)
legend('4 PPC','6 PPC','8 PPC')
title('RMS error with volume updating','Fontsize',12)

% Number of elements
el = [5 10 15 20 25 30 35 ];

% Error at time t = 5, delta t = 0.001
er = [0.4906 0.1294 0.0578 0.0326 0.0208 0.0145 0.0106];

% Error at time t = 5 , delta t = t_cr * CFL 
er_half = [0.4757 0.1279 0.0570 0.0324 0.0207 0.0143 0.0106];

% Error at time t = 0.5, delta t = 0.001
er2 = [0.0644 0.0167 0.0074  0.0042 0.0027 0.0019 0.0014 ];

% Error at time t = 0.5, delta t = t_cr * CFL
er_half2 = [0.0613 0.0160 0.0070 0.0041 0.0026 0.0018 0.0013]


er_onethird = [0.4217 0.1123 0.0501 0.0284 0.0182 0.0126 0.0093];


figure(7)
plot(log2(el),log2(er))
xlabel('Log(Number of elements)')
ylabel('Log(Error)')
 
figure(8)
plot(el,er)
xlabel('Number of elements')
ylabel('Error')

%% File in which the error of the MPM 1D solution of the vibrating string solution is stored.

% Number of elements
el = [5 10 15 20 25 30 35 ];

% 1 particle per elements Error at time t = 0.5, delta t = 0.001. Length of simulation 5 sec.
er = [0.0029 7.3745e-04 3.2869e-04 1.8504e-04 1.1845e-04 8.2247e-05    6.0410e-05];

% 2 particles per element Error at time t = 0.5, delta t = 0.001. Length of simulation 5 sec.  
er2 = [0.0026    6.5349e-04 2.9127e-04  1.6398e-04 1.0497e-04 7.2887e-05 5.3536e-05]

% 4 particles per element Error at time t = 0.5, delta t = 0.001. Length ofsimulation 5 sec
er4 = [0.0023 5.8996e-04 2.6295e-04  1.4803e-04 9.4760e-05 6.5802e-05 4.8334e-05]

% 10 particles per element Error at time t = 0.5, delta t = 0.001. Length ofsimulation 5 sec
er10 = [0.0020 5.1800e-04 2.2804e-04 1.3226e-04  8.3796e-05 5.8465e-05 4.3509e-05  ]

%figure(1)
%plot(log2(el),log2(er))
%xlabel('Log(Number of elements)')
%ylabel('Log(Error)')
 
figure(9)
plot(el,er)

hold on
% Plot of the error for different number of particles per element and
% different number of elements 
plot(el,er2,'r')
plot(el,er4,'g')
plot(el,er10,'m')
xlabel('Number of elements')
ylabel('Error')
legend('1 part/el','2 part/el','4 part/el','10 part/el')

%% File in which the error of the MPM 1D solution of the oedometer test solution is stored.
%% Input 
% density = 1;
% Youngs_modulus = 50000;
% gravitational_acceleration = -9.8; 
% load = 0;
% height = 25;

% Number of elements
el = [10 15 20 25 30 35 40 45 50 55 60 65 75 85];

% 4 particles per element Error at time t = 0.5, delta t = 0.001. Length of simulation 2.5 sec. With 40 elements grid crossing occured!
vol_er4 =  [1.4038e-03 5.3587e-04 3.3903e-04 2.5076e-04 1.7161e-04 1.3086e-04 1.1077e-04 6.1505e-04 1.7285e-03 1.8934e-03 2.7263e-03 2.9518e-03 4.7722e-03 6.9001e-03 ];          

% 6 particles per element Error at time t = 0.5, delta t = 0.001. Length of simulation 2.5 sec. With 30 elements grid crossing occured! 
vol_er6 =  [1.3517e-03 5.1419e-04 3.3277e-04 2.4293e-04 1.8990e-04 9.2701e-04 1.8257e-03 1.4905e-03 2.4439e-03 2.9544e-03 1.6745e-03 2.705e-03 2.7579e-03 5.004e-03];

% 8 particles per element Error at time t = 0.5, delta t = 0.001. Length of simulation 2.5 sec. With 25 elements grid crossing occured!
vol_er8 =  [1.3214e-03 5.0285e-04 3.2917e-04 3.6547e-04 7.0347e-04 1.1464e-03 1.7413e-03 1.2723e-03 1.2038e-03 2.2279e-03 1.6085e-03 2.1298e-03 1.7715e-03 3.6589e-03];  

% Particle boundary crossing starts at 80 elements 

figure(10)
semilogy(el,vol_er4,'-*r','LineWidth',2)
hold on
semilogy(el,vol_er6,'-*','LineWidth',2)
semilogy(el,vol_er8,'-*m','LineWidth',2)
xlabel('Number of elements','Fontsize',12)
ylabel('RMS Error','Fontsize',12)
legend('4 PPC','6 PPC','8 PPC','Location','northwest')
title('RMS error with increasing PPC','Fontsize',12)

%% File in which the error of the MPM 1D solution of the vibrating string solution is stored.
%% Input
% density = 50; 
% Youngs_modulus = 1E5;
% gravitational_acceleration = -9.81; 
% load = -10;
% height = 1;
% number_particles_per_element = 4;
% t_step = 1E-3; 


% Number of elements
el = [5 10 15 20 25 30 35 36 37];

% RME error at t = 0.05 sec
Error = [ 5.0885e-05 4.3186e-05  4.0865e-05 3.9945e-05 3.9731e-05  3.9800e-05 3.9874e-05  3.9882e-05 3.9948e-05];

% RME error at t = 0.5 sec
Error2 = [1.6033e-04 3.7690e-05 1.0688e-04 1.2303e-04 1.2971e-04  1.3454e-04 1.3841e-04  1.3905e-04  2.9678e-04 ];

figure(11)
semilogy(el, Error,'-*','LineWidth',2)
axis([0 40 10E-7 10E-4])
xlabel('Number of elements','Fontsize',12)
ylabel('RMS Error','Fontsize',12)
legend('4 PPC')
title('RME error at t = 0.05 s','Fontsize',12)
% 38 elements -> particle cross boundary -> Instability

figure(12)
semilogy(el, Error2,'-*','LineWidth',2)
xlabel('Number of elements','Fontsize',12)
ylabel('RMS Error','Fontsize',12)
legend('4 PPC')
title('RME error at t = 0.5 s','Fontsize',12)

%% File in which the error of the MPM 1D solution of the vibrating string solution is stored.
% Input
% density = 1;
% Youngs_modulus = 100;
% gravitational_acceleration = 0; 
% load = 0;
% height = 25; 

% Number of elements
el = [5 10 15 20 25 30 35 ];

% 1 particle per elements Error at time t = 0.5, delta t = 0.001. Length of simulation 5 sec.
er = [0.0029 7.3745e-04 3.2869e-04 1.8504e-04 1.1845e-04 8.2247e-05    6.0410e-05];

% 2 particles per element Error at time t = 0.5, delta t = 0.001. Length of simulation 5 sec.  
er2 = [0.0026    6.5349e-04 2.9127e-04  1.6398e-04 1.0497e-04 7.2887e-05 5.3536e-05]

% 4 particles per element Error at time t = 0.5, delta t = 0.001. Length ofsimulation 5 sec
er4 = [0.0023 5.8996e-04 2.6295e-04  1.4803e-04 9.4760e-05 6.5802e-05 4.8334e-05]

% 10 particles per element Error at time t = 0.5, delta t = 0.001. Length of simulation 5 sec
er10 = [0.0020 5.1800e-04 2.2804e-04 1.3226e-04  8.3796e-05 5.8465e-05 4.3509e-05  ]


figure(13)
semilogy(el,er4,'-*','LineWidth',2)
xlabel('Number of elements','Fontsize',12)
ylabel('RMS Error','Fontsize',12)
legend('4 PPC')

figure(14)
semilogy(el,er,'-*','LineWidth',2)
hold on
semilogy(el,er2,'-*r','LineWidth',2)
semilogy(el,er4,'-*g','LineWidth',2)
semilogy(el,er10,'-*m','LineWidth',2)
xlabel('Number of elements','Fontsize',12)
ylabel('RMS Error','Fontsize',12)
legend('1 PPC','2 PPC','4 PPC','10 PPC')

%% Test problem steffen: spatial convergence

% 3.5 PPC ,simulation 4 sec, hyper elasticity , CFL = 0.5, reset 
% error at t = 1 s. explicit euler forward.
elementen = [10 20 40 80 160 320]
RME_error = [ 3.0945e-04   3.8832e-04  4.4959e-04   4.8569e-04    4.9891e-04    5.0547e-04 ]

figure(15)
loglog(elementen,RME_error,'-*','LineWidth',2)
xlabel('Number of elements','Fontsize',12)
ylabel('RMS Error','Fontsize',12)
axis([10 10e2 10e-8 10e-3])
title('Spatial convergence MPM','Fontsize',12)
legend('MPM')
