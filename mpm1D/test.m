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


figure(1)
plot(log2(el),log2(er))
xlabel('Log(Number of elements)')
ylabel('Log(Error)')
 
figure(2)
plot(el,er)
xlabel('Number of elements')
ylabel('Error')

