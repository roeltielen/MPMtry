function [x,y] = exact_solution_dynamic_bar(rho, tau, l, x0, t, number_time_steps)
% Input
% t = 0:0.01:1.5;
% rho = 10^2; 
% tau = 1; 
% length = 1;
% Initial position x0

% Parameters
alpha = l*tau/(rho*pi);
omega = pi/l;
x = zeros(1,number_time_steps);
y = zeros(1,number_time_steps);

%% Analytical solution for displacement
for i = 1:1:number_time_steps
    if (t(i) < l - x0)
        x(1,i) = 0;
        y(1,i) = 0;
    end
        if (l - x0 <= t(i)) && (t(i) < l + x0)
        x(1,i) = alpha*(1 + cos(omega*(t(i)+x0)));
        y(1,i) = -alpha*omega*sin(omega*(t(i)+x0));
        end
            if (l + x0 <= t(i)) &&(t(i) < 3*l - x0)
            x(1,i) = alpha*(cos(omega*(t(i)+x0))-cos(omega*(t(i)-x0)));
            y(1,i) = -alpha*omega*(sin(omega*(t(i)- x0))+ sin(omega*(t(i)+ x0)));
            end
                if (3*l - x0 <= t(i)) && (t(i) < 3*l + x0)
                x(1,i) = alpha*(-1-cos(omega*(t(i)-x0)));
                y(1,i) = alpha*omega*sin(omega*(t(i)-x0));
               
                end                   
end
x;
end