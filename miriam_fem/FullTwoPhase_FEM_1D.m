%% Miriam Mieremet, 10-March-2015
% Full 2-phase
% 
% rw * w_t = p_z - rw * g - (n * rw * g / k) * (w - v)
% (1 - n) * rs * v_t = - n * rw * w_t + sigma_z + p_z - rsat * g 
% p_t = (Kw / n) * [(1 - n) * v_z + n * w_z]
% sigma_t = E * v_z
% u_t = v

% 1D TL-FEM Euler-Cromer

clear all
clc

%% Get input
% Read all input from excelsheet

nod = xlsread('FullTwoPhase_Input1D.xlsx','Nodes');
n_nod = length(nod);
H = max(nod);

elm = xlsread('FullTwoPhase_Input1D.xlsx','Elements');
n_elm = length(elm);

mat = xlsread('FullTwoPhase_Input1D.xlsx','Material');
rs = mat(1);    % Density Solid
Ec = mat(2);    % Young's Modulus
nu = mat(3);    % Poisson Ratio
n = mat(4);     % Porosity
k = mat(5);     % Hydraulic conductivity
rw = mat(6);    % Density Fluid
Kw = mat(7);    % Bulk Modulus

ic = xlsread('FullTwoPhase_Input1D.xlsx','InitialCond');
u = ic(:,1);
v = ic(:,2);
w = ic(:,3);
p = ic(1:n_elm,5);
s = ic(1:n_elm,6);

fix = xlsread('FullTwoPhase_Input1D.xlsx','Fixities');
n_fix = length(fix(:,1));

load = xlsread('FullTwoPhase_Input1D.xlsx','Load');
n_belm = length(load(:,1));

time = xlsread('FullTwoPhase_Input1D.xlsx','TimeInt');
dt = time(1);   % Size of time step
m = time(2);    % Number of time steps

% Close all Excel processes in task manager
% system('taskkill /F /IM EXCEL.EXE'); 

%% Define other variables

g = 10;                                        % Gravitational acceleration
rsat = (1 - n) * rs + n * rw;                  % Density of mixture
gammaw = rw * g;                               % Unit weight fluid
gamma = rsat * g;                              % Unit weight mixture

%% Define shape function for elements
% First order interpolation functions
% Single Gauss Point in normalized linear element

N = 0.5;
LN = [-1,1];

%% Define shape function for boundary elements
% First order interpolation functions
% Single Gauss Point in normalized point element

N_b = 1;

%% Define lumped mass matrix of fluid and solid (time-independent)

Mw = zeros(n_nod,1);
Ms = zeros(n_nod,1);

% for-loop over elements
for i=1:n_elm
    getnod = elm(i,:);
    J = nod(getnod(2)) - nod(getnod(1));
    
    Mw(getnod(1)) = Mw(getnod(1)) + N * rw * J;
    Mw(getnod(2)) = Mw(getnod(2)) + N * rw * J;
    
    Ms(getnod(1)) = Ms(getnod(1)) + N * (1 - n) * rs * J;
    Ms(getnod(2)) = Ms(getnod(2)) + N * (1 - n) * rs * J;
end

M_w = n * Mw;

%% Define traction force vector for fluid and mixture (time-independent)

F_trac_w = zeros(n_nod,1);
F_trac = zeros(n_nod,1);

% for-loop over boundary elements
for i=1:n_belm
    getload = load(i,:);
        
    if(getload(4)==1)
        F_trac_w(getload(1)) = F_trac_w(getload(1)) + N_b * getload(5);
    end
    
    if(getload(2)==1)
        F_trac(getload(1)) = F_trac(getload(1)) + N_b * getload(3);
    end
end

%% Define gravity force vector for fluid and mixture (time-independent)

% F_grav_w = zeros(n_nod,1);
% F_grav = zeros(n_nod,1);
% 
% % for-loop over elements
% for i=1:n_elm
%     getnod = elm(i,:);
%     J = nod(getnod(2)) - nod(getnod(1));
%     
%     F_grav_w(getnod(1)) = F_grav_w(getnod(1)) - N * gammaw * J;
%     F_grav_w(getnod(2)) = F_grav_w(getnod(2)) - N * gammaw * J;
%     
%     F_grav(getnod(1)) = F_grav(getnod(1)) - N * gamma * J;
%     F_grav(getnod(2)) = F_grav(getnod(2)) - N * gamma * J;
% end

%% Define lumped drag matrix (time-independent)

Q = zeros(n_nod,1);

% for-loop over elements
for i=1:n_elm
    getnod = elm(i,:);
    J = nod(getnod(2)) - nod(getnod(1));
    Q(getnod(1)) = Q(getnod(1)) + N * n * gammaw * J / k;
    Q(getnod(2)) = Q(getnod(2)) + N * n * gammaw * J / k;
end

%% TIME INTEGRATION

saveu = zeros(n_nod,m+1);
saveu(:,1) = nod + u;
%savev = zeros(n_nod,m+1);
%savew = zeros(n_nod,m+1);
savep = zeros(n_elm,m+1);
savep(:,1) = p;
saves = zeros(n_elm,m+1);
saves(:,1) = s;


for t = 1:m
    %% Define internal force vector for effective stress and pore pressure (time-dependent)

    F_int_s = zeros(n_nod,1);
    F_int_p = zeros(n_nod,1);
    
    % for-loop over elements
    for i=1:n_elm
        getnod = elm(i,:);
        F_int_s(getnod(1)) = F_int_s(getnod(1)) + LN(1)*s(i);   % effective stress
        F_int_s(getnod(2)) = F_int_s(getnod(2)) + LN(2)*s(i);
        
        F_int_p(getnod(1)) = F_int_p(getnod(1)) + LN(1)*p(i);   % pore pressure
        F_int_p(getnod(2)) = F_int_p(getnod(2)) + LN(2)*p(i);
    end
       
    %% Update velocity, effective stress and displacement
    F_dr = Q .* (w - v);
    
    % update velocity fluid (explicit)
    a = Mw .\ (- F_int_p + F_trac_w - Q .* (w - v)); % + F_grav_w zonder gravitatie
    
    % (fixities)
    for i = 1:n_fix
        if(fix(i,2)==1)
            a(fix(i,1)) = 0;
        end
    end
    
    w = w + a .* dt;
    
    % update velocity solid (explicit)
    v = v + Ms .\ ( - M_w .* a + F_trac - F_int_s - F_int_p ) .* dt; 
    % + F_grav zonder gravitatie 
 
    % (fixities)
    for i = 1:n_fix
        if(fix(i,2)==1)
            v(fix(i,1)) = 0;
        end
    end
    
    % update pore pressure (implicit)
    dpdt = zeros(n_elm,1);
    for i=1:n_elm
        getnod = elm(i,:);
        J = nod(getnod(2)) - nod(getnod(1));
        dpdt(i) = (Kw / n) * ((1 - n) * (LN(1) * v(getnod(1)) + LN(2) * v(getnod(2))) + n * (LN(1) * w(getnod(1)) + LN(2) * w(getnod(2)))) / J;
    end
    
    p = p + dt * dpdt;  
    
    % update stress (implicit)
    
    % (calculate strain rate)
    dedt = zeros(n_elm,1);
    for i=1:n_elm
        getnod = elm(i,:);
        J = nod(getnod(2)) - nod(getnod(1));
        dedt(i) = (LN(1) * v(getnod(1)) + LN(2) * v(getnod(2)))/J;
    end
    
    s = s + dt * Ec * dedt; 
    
    % update displacement (implicit)
    u = u + dt .* v;
    
    saveu(:,t+1) = nod + u;
    %savev(:,t+1) = v;
    %savew(:,t+1) = w;
    savep(:,t+1) = p;
    saves(:,t+1) = s;
end

hold off

% plot numerical solution for x = 1.5
t = 0 : dt : m*dt;
plot(t,saves(601,:)/(load(5)),'-','Color',...
    [0.301960796117783 0.745098054409027 0.933333337306976],...
    'LineWidth',1)
hold on

% plot analytical solution for x = 1.5

t = 0 : dt/4 : m*dt;



c = sqrt((Ec)/rsat);
cv = k/rw/g/(1/(Ec)+1/(Kw/n));
df = rw/rsat;
ds = rs/rsat;
b = (Ec)/(Kw/n);

T = 30*H/sqrt((Ec+Kw/n)/rsat);
x = 0.5;
K = 1001;
AS = 17/100 + 0 * x;
for j = 1:2:K
    j
    omega = 2*pi*j/T;
    a = n*g/k/omega;

    A = n;
    B = -n*(1-n)*ds - (1-n)^2*df - df*b + complex(0,1)*a*df*(1+b);
    C = (1-n)*ds*df*b - complex(0,1)*a*df*b;
    
    gam1 = -sqrt((-B + sqrt(B^2-4*A*C))/2/A);
    gam2 = -sqrt((-B - sqrt(B^2-4*A*C))/2/A);
    
    Ap = (df*b-gam2^2)*(1-df-gam1^2)/(gam1^2-gam2^2)/(1-df-b*df)*2/pi/j;
    Bp = -(df*b-gam1^2)*(1-df-gam2^2)/(gam1^2-gam2^2)/(1-df-b*df)*2/pi/j;
    
    As = -(df*b-gam1^2)/(1-df-gam1^2) * Ap;
    Bs = -(df*b-gam2^2)/(1-df-gam2^2) * Bp;
    
    AS = AS + imag(As * exp(-(omega/c)*(imag(gam1)-complex(0,1)*real(gam1))*x)*exp(complex(0,1)*omega*t) + Bs * exp(-(omega/c)*(imag(gam2)-complex(0,1)*real(gam2))*x)*(exp(complex(0,1)*omega*t))...
    + As * exp(-(omega/c)*(imag(gam1)-complex(0,1)*real(gam1))*(2*H-x))*exp(complex(0,1)*omega*t) + Bs * exp(-(omega/c)*(imag(gam2)-complex(0,1)*real(gam2))*(2*H-x))*(exp(complex(0,1)*omega*t))...
    - As * exp(-(omega/c)*(imag(gam1)-complex(0,1)*real(gam1))*(x+2*H))*exp(complex(0,1)*omega*t) - Bs * exp(-(omega/c)*(imag(gam2)-complex(0,1)*real(gam2))*(x+2*H))*(exp(complex(0,1)*omega*t)));    
end

plot(t,AS,'-k','LineWidth',1)

c_1 = sqrt((Ec+Kw/n)/rsat)
d_1 = 0.50/c_1
e_1 = 3.50/c_1

c_2 = sqrt((Ec*n)/((1-n)*Kw+n*Ec))*sqrt((Kw)/rw)
d_2 = 0.50/c_2
e_2 = 3.50/c_2

% plottools 
xlim([0 m*dt])
ylim([0.0 2.0])
xlabel('time (s)','FontSize',28)
ylabel('normalized effective stress (-)','FontSize',28)
set(gca,'YTickLabel',{'0','0.5','1.0','1.5','2.0'},...
    'YTick',[0 0.5 1 1.5 2],...
    'XTickLabel',{'0','0.001','0.002','0.003'},...
    'XTick',[0 0.001 0.002 0.003],...
    'FontSize',24)
l = legend('numerical solution','analytical solution',...
    'Location','NorthWest');
set(l,'FontSize',28);






