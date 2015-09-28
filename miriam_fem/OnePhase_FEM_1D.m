%% Miriam Mieremet, 01-July-2014
% 1-phase
% 
% rs * v_t = sigma_z + rs * g
% sigma_t = E * v_z
% u_t = v
% 
% 1D TL-FEM Euler-Cromer

clear all
clc

%% Get input
% Read all input from excelsheet

nod = xlsread('OnePhase_Input1D.xlsx','Nodes');
n_nod = length(nod);
H = max(nod);

elm = xlsread('OnePhase_Input1D.xlsx','Elements');
n_elm = length(elm);

mat = xlsread('OnePhase_Input1D.xlsx','Material');
rho = mat(1);    % Density Solid
E = mat(2);      % Young's Modulus
nu = mat(3);     % Poisson Ratio

ic = xlsread('OnePhase_Input1D.xlsx','InitialCond');
u = ic(:,1);
v = ic(:,2);
s = ic(1:n_elm,5);

fix = xlsread('OnePhase_Input1D.xlsx','Fixities');
n_fix = length(fix(:,1));

load = xlsread('OnePhase_Input1D.xlsx','Load');
n_load = length(load(:,1));
p = load(1,2);

time = xlsread('OnePhase_Input1D.xlsx','TimeInt');
dt = time(1);   % Size of time step
m = time(2);    % Number of time steps
a = time(3);    % Damping factor

% Close all Excel processes in task manager
% system('taskkill /F /IM EXCEL.EXE');  

%% Define other variables

g = 10;        % Gravitational acceleration

%% Define shape function for elements
% First order interpolation functions
% Single Gauss Point in normalized linear element

N = 0.5;
LN = [-1,1];

%% Define shape function for boundary elements
% First order interpolation functions
% Single Gauss Point in normalized point element

N_b = 1;

%% Define lumped-mass matrix (time-independent)

M = zeros(n_nod,1);

% for-loop over elements
for i=1:n_elm
    getnod = elm(i,:);
    J = nod(getnod(2)) - nod(getnod(1));
    M(getnod(1)) = M(getnod(1)) + N * rho * J;
    M(getnod(2)) = M(getnod(2)) + N * rho * J;
end

%% Define traction force vector (time-independent)

F_trac = zeros(n_nod,1);

% for-loop over boundary elements
for i=1:n_load
    getload = load(i,:);
    F_trac(getload(1)) = F_trac(getload(1)) + N_b * getload(2);
end

%% Define gravity force vector (time-independent)

F_grav = - M .* g;


%% TIME INTEGRATION

saveu = zeros(n_nod,m+1);
saveu(:,1) = nod;
savev = zeros(n_nod,m+1);
saves = zeros(n_elm,m+1);

for t = 1:m
    %% Define internal force vector (time-dependent)

    F_int = zeros(n_nod,1);

    % for-loop over elements
    for i=1:n_elm
        getnod = elm(i,:);
        F_int(getnod(1)) = F_int(getnod(1)) + LN(1)*s(i);
        F_int(getnod(2)) = F_int(getnod(2)) + LN(2)*s(i);
    end
    
    %% Update velocity, stress and displacement
    
    % update velocity
    v = v + M .\ (F_trac - F_int + F_grav) .* dt;
    
    % fixities
    for i = 1:n_fix
        if(fix(i,2)==1)
            v(fix(i,1)) = 0;
        end
    end
    
    % calculate strain rate
    de = zeros(n_elm,1);
    for i=1:n_elm
        getnod = elm(i,:);
        J = nod(getnod(2)) - nod(getnod(1));
        de(i) = (LN(1) * v(getnod(1)) + LN(2) * v(getnod(2)))/J;
    end
    
    % update stress
    s = s + dt * E * de;
    
    % update displacement
    u = u + dt * v;    
    
    saveu(:,t+1) = nod + u;
    savev(:,t+1) = v;
    saves(:,t+1) = s;
end

hold off

figure1 = figure;
set(figure1,'units','normalized','outerposition',[0 0 1 1])
axes1 = axes('Parent',figure1,...
    'Position',[0.13 0.11 0.524404761904762 0.815],...
    'FontSize',20,...
    'FontName','Times New Roman');
xlim(axes1,[0 0.5]);
%ylim(axes1,[0 1]);
box(axes1,'on');
hold(axes1,'all');

% plot numerical solution for x = {H, 3H/4, H/2, 0}
t = 0 : dt : m*dt;
plot(t,saveu(n_nod,:),'o','MarkerSize',4,'Color',[1 0 0],'DisplayName','$x_0 = 1.00$')
hold on
plot(t,saveu((n_nod-1)*3/4+1,:),'o','MarkerSize',4,'Color',[0 0 1],'DisplayName','$x_0 = 0.75$')
plot(t,saveu((n_nod-1)/2+1,:),'o','MarkerSize',4,'Color',[0 1 0],'DisplayName','$x_0 = 0.50$')
plot(t,saveu((n_nod-1)/4+1,:),'o','MarkerSize',4,'Color',[1 0 1],'DisplayName','$x_0 = 0.25$')

% plot exact solution for x = {H, 3H/4, H/2, 0}
t =  0 : dt/100 : m*dt;
UH44 = (1/2)*rho*g*H^2/E+(p-rho*g*H)*H/E-(8*pi*p-...
    16*rho*g*H)*H*cos((1/2)*sqrt(E/rho)*pi*t/H)/(pi^3*E)+...
    (1/27)*(-24*pi*p-16*rho*g*H)*H*cos((3/2)*sqrt(E/rho)*pi*t/H)/(pi^3*E)-...
    (1/125)*(40*pi*p-16*rho*g*H)*H*cos((5/2)*sqrt(E/rho)*pi*t/H)/(pi^3*E)+...
    (1/343)*(-56*pi*p-16*rho*g*H)*H*cos((7/2)*sqrt(E/rho)*pi*t/H)/(pi^3*E);
UH34 = (9/32)*rho*g*H^2/E+(3/4)*(p-rho*g*H)*H/E-...
    (8*pi*p-16*rho*g*H)*H*cos((1/2)*sqrt(E/rho)*pi*t/H)*sin((3/8)*pi)/(pi^3*E)+...
    (1/27)*(-24*pi*p-16*rho*g*H)*H*cos((3/2)*sqrt(E/rho)*pi*t/H)*sin((1/8)*pi)/(pi^3*E)+...
    (1/125)*(40*pi*p-16*rho*g*H)*H*cos((5/2)*sqrt(E/rho)*pi*t/H)*sin((1/8)*pi)/(pi^3*E)-...
    (1/343)*(-56*pi*p-16*rho*g*H)*H*cos((7/2)*sqrt(E/rho)*pi*t/H)*sin((3/8)*pi)/(pi^3*E);
UH24 = (1/8)*rho*g*H^2/E+(1/2)*(p-rho*g*H)*H/E-...
    (1/2)*(8*pi*p-16*rho*g*H)*H*cos((1/2)*sqrt(E/rho)*pi*t/H)*sqrt(2)/(pi^3*E)-...
    (1/54)*(-24*pi*p-16*rho*g*H)*H*cos((3/2)*sqrt(E/rho)*pi*t/H)*sqrt(2)/(pi^3*E)+...
    (1/250)*(40*pi*p-16*rho*g*H)*H*cos((5/2)*sqrt(E/rho)*pi*t/H)*sqrt(2)/(pi^3*E)+...
    (1/686)*(-56*pi*p-16*rho*g*H)*H*cos((7/2)*sqrt(E/rho)*pi*t/H)*sqrt(2)/(pi^3*E);
UH14 = (1/32)*rho*g*H^2/E+(1/4)*(p-rho*g*H)*H/E-...
    (8*pi*p-16*rho*g*H)*H*cos((1/2)*sqrt(E/rho)*pi*t/H)*sin((1/8)*pi)/(pi^3*E)-...
    (1/27)*(-24*pi*p-16*rho*g*H)*H*cos((3/2)*sqrt(E/rho)*pi*t/H)*sin((3/8)*pi)/(pi^3*E)-...
    (1/125)*(40*pi*p-16*rho*g*H)*H*cos((5/2)*sqrt(E/rho)*pi*t/H)*sin((3/8)*pi)/(pi^3*E)-...
    (1/343)*(-56*pi*p-16*rho*g*H)*H*cos((7/2)*sqrt(E/rho)*pi*t/H)*sin((1/8)*pi)/(pi^3*E);
plot(t,4*H/4+UH44,'k','HandleVisibility','off')
plot(t,3*H/4+UH34,'k','HandleVisibility','off')
plot(t,2*H/4+UH24,'k','HandleVisibility','off')
plot(t,1*H/4+UH14,'k','HandleVisibility','off')

% Create xlabel
xlabel('$t$ ($s$)','Interpreter','latex','FontSize',24,...
    'FontName','Times New Roman');

% Create ylabel
ylabel('$x$ ($m$)','Interpreter','latex','FontSize',24,...
    'FontName','Times New Roman');

% Create text
text('Parent',axes1,'Interpreter','latex',...
    'String',{['$\begin{array}{lcrl} \rho &=&  ',num2str(rho),'& kg/m^3\\ E &=& ',...
    num2str(E),'& Pa\\ g &=& ', num2str(g),'& m/s^3\quad \\ \\ H &=& ',...
    num2str(H),'& m\\ p_0 &=& ', num2str(p),'& Pa \\ \\ \Delta t &=&',...
    num2str(dt),'& s\\ \Delta x_0 &=&',num2str(H/n_elm),'& m\end{array}$']},...
    'Position',[0.512365002365001 0.226085525138741 0],...
    'FontSize',24,...
    'EdgeColor',[0 0 0],...
    'BackgroundColor',[1 1 1]);

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Interpreter','latex',...
    'Position',[0.735714285714286 0.678932176368448 0.20952380952381 0.245955716008234],...
    'FontSize',24);





% hold off
% for j = 1: m
%     subplot(1,3,1);
%     plot(0*saveu(:,j)+1,saveu(:,j),'-x')
%     axis([0 2 0 1])
%     
%     subplot(1,3,2);
%     plot(savev(:,j),saveu(:,j))
%     axis([min(min(savev)) max(max(savev)) 0 1])
%     
%     subplot(1,3,3);
%     plot(saves(:,j),(saveu(1:n_nod-1,j)+saveu(2:n_nod,j))/2)
%     axis([min(min(saves)) max(max(saves)) 0 1])
%     
%     frame(j) = getframe(gcf);
% end
% 
% movie(gcf,frame,1,5)
