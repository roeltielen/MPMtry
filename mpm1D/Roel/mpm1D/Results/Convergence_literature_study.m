%% Spatial oedometer problem (adjusted parameters)
% Deformation = 1 and reset = 1 !
% Used parameters:
% Density = 1
% E = 5E4
% g = -9.81 
% Length = 25
% Delta t = 1E-3
close all

% Length of simulation 2.5 , RMS error at 0.5 s
el = [ 4 8 12 16 24 32 48 64 96 128];

er4 = [0.0041 0.0019  7.7723e-04      3.8530e-04  3.5351e-04   3.4308e-04 0.0013   0.0029  0.0062  0.0164];

er6 = [0.0040 0.0018  7.4299e-04      3.7734e-04  3.5202e-04   5.8636e-04 0.0028   0.0025  0.0042  0.0084];

er8 = [0.0040 0.0018  7.2461e-04      3.7349e-04  4.3493e-04   0.0011  0.0022     0.0027  0.0042  0.0063];

figure(2)
set(gcf, 'PaperPosition', [0 0 7 7]);
set(gcf, 'PaperSize', [6 4.5]);
semilogy(el,er4,'-*r','LineWidth',2)
hold on
semilogy(el,er6,'-*g','LineWidth',2)
semilogy(el,er8,'-*','LineWidth',2)
xlabel('Number of elements','Fontsize',12)
ylabel('RMS Error','Fontsize',12)
legend('4 PPC','6 PPC','8 PPC')
title('Spatial convergence linear MPM','Fontsize',12)
axis([0 128 10e-5 10e-2])


%% Spatial convergence vibrating string (both ends fixed)
% Deformation = 1 and reset = 1 !
% Used parameters:
% Density = 1
% E = 100
% g = 0 
% v_0 = 0.1
% Length = 25
% Delta t = 1E-3

% Length of simulation 2.5 , RMS error at 0.5 s
el = [ 4 8 16 32 64];

er4 = [0.0036 9.1906e-04 2.3116e-04 5.7897e-05 1.4756e-05];

er6 = [0.0034 8.7821e-04 2.2087e-04 5.5329e-05 1.4134e-05];

er8 = [0.0033 8.5591e-04 2.1525e-04 5.3926e-05 1.3831e-05];

figure(1)
set(gcf, 'PaperPosition', [0 0 7 7]);
set(gcf, 'PaperSize', [6 4.5]);
semilogy(el,er4,'-*r','LineWidth',2)
hold on
semilogy(el,er6,'-*g','LineWidth',2)
semilogy(el,er8,'-*','LineWidth',1.5)
xlabel('Number of elements','Fontsize',12)
ylabel('RMS Error','Fontsize',12)
legend('4 PPC','6 PPC','8 PPC')
title('Spatial convergence linear MPM','Fontsize',12)
axis([0 64 10e-6 10e-3])