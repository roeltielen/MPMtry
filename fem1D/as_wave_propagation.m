function[t,AS] = as_wave_propagation()
% Aangepaste code van Miriam
% plot analytical solution for x = 1.5
Ec = 2.5E9; %5E9;
rw = 1E3;
rs = 2.6E3;
n = 0.4;
rsat = (1-n)*rs + n*rw;
k = 1E-5;
g = 10;
Kw = 2E9;
H = 2;

dt = 1E-5;
m = 400;
t = 0 : dt/4 : m*dt;

 
 
c = sqrt((Ec)/rsat);
cv = k/rw/g/(1/(Ec)+1/(Kw/n));
df = rw/rsat;
ds = rs/rsat;
b = (Ec)/(Kw/n);

T = 30*H/sqrt((Ec+Kw/n)/rsat);
x = 0.5;
K = 10001;
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

figure(1);
plot(t,AS,'-k','LineWidth',1)
