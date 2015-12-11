[node_4,pos_4,t] = oedometer_1phase(4);
[node_8,pos_8] = oedometer_1phase(8);
[node_16,pos_16] = oedometer_1phase(16);

figure(15)
plot(t,pos_4,'r','LineWidth',2)
hold on
plot(t,pos_8,'m','LineWidth',2)
hold on
plot(t,pos_16,'b','LineWidth',2)
% plot(t,position_exact_nodes(node_number+1,:),'--r','LineWidth',2)
xlabel('time [s]','Fontsize',12)
ylabel('position [m]', 'Fontsize',12)
%title(sprintf('Position of node %d',node_4),'FontSize', 12)
legend('h=1/4','h=1/8','h=1/16')
hold on
set(gcf, 'PaperUnits','normalized')
set(gcf, 'PaperPosition', [0 0 1 1]);
%set(gcf, 'PaperSize', [6 4.7]);