function [position_exact] = exact_solution(pos_p_glob,n,t_step,v_0,Youngs_modulus,...
    density,height)
%% Obtain the exact solution for the particles(!)
number_particles = size(pos_p_glob,1);
position_exact = zeros(number_particles, 1);
disp_exact = zeros(number_particles, 1);
%vel_exact = zeros(number_particles, 1);


w1 = pi*sqrt(Youngs_modulus/density)/(2*height);
b1 = pi/(2*height);
for p = 1:number_particles
    %vel_exact(p,1) = v_0*cos(w1*t_step*n)*sin(b1*pos_p_glob(p));
    disp_exact(p,1) = v_0/w1*sin(w1*t_step*n)*sin(b1*pos_p_glob(p));
end


position_exact(:,1) = disp_exact(:,1) + pos_p_glob;


end