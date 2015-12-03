%% Representation of the results regarding numerical integration
% f(x) = x^2 on [0,1] 
% 1 Element, 1-2-3 PPC

node_1_1 = 1/8;
node_1_2 = 1/8; 

node_2_1 = 1/9;
node_2_2 = 1/6;

node_3_1 = 5/48;
node_3_2 = 9/48;

node_4_1 = 1/10;
node_4_2 = 1/5;

node_exact_1 = 1/12;
node_exact_2 = 1/4;

plot(node_1_1,node_1_2,'*r')
hold on
plot(node_2_1,node_2_2,'*k')
plot(node_3_1,node_3_2,'*m')
plot(node_4_1,node_4_2,'*g')
plot(node_exact_1,node_exact_2,'*')
xlabel('Projection on node 1')
ylabel('Projection on node 2')
title('Illustration quality MPM integration')
legend('1 PPC', '2 PPC', '3 PPC', '4 PPC','Exact')
axis([0 0.25 0 0.4])
