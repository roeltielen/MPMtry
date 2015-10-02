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
 
figure(2)
plot(el,er)

hold on
plot(el,er2,'r')
plot(el,er4,'g')
plot(el,er10,'m')
xlabel('Number of elements')
ylabel('Error')
legend('1 part/el','2 part/el','4 part/el','10 part/el')