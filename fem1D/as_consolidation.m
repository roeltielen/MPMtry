function[x,AS] = as_consolidation()
% Aangepaste code van Miriam

Ec = 2.5E9; %1E7; %5E9;
rw = 1E3;
rs = 2600; %(1/0.6)*1.6E3;
n = 0.4;
%rsat = (1-n)*rs + n*rw;
k = 1E-5;
g = 10;
Kw = 2E9; %3E8;
H = 1;
J = 1/800;

% plot analytical solution
cv = k / (rw * g * (1/Ec + n/Kw)); 
T = [0.02,0.05,0.1,0.2,0.5,1];
x = (0:J/2:H+0.00000001)';
K = 50001;
AS = 0*x*T;
for j = 1:K
    for i = 1:length(T)
    AS(:,i) = AS(:,i) + 4/pi * (-1)^(j-1)/(2*j-1) * cos((2*j-1) * pi/2 * (x) / H) * exp(-(2*j-1)^2 * pi^2/4 * T(i));
    end
end
figure(5)
for i = 1:length(T)
    plot(AS(:,i),x,'-k','LineWidth',1)
    hold on
end
