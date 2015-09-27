function[M_loc, F_grav_loc, F_int_loc, B_loc] = gauss_rule(density, E, h, g)
%..Note: all computations take place in the local domain
 
%..Prescribe local positions and weights of Gauss points
pos_loc = 0.5;
weight  = 1;

%..Define basis functions vector 
N_loc = [1-pos_loc, pos_loc];

%..Define strain-displacement matrix 
B_loc = [-1/h, 1/h];

%..Compute mass matrix  
M_loc = weight*N_loc'*density*N_loc*h;

%..Compute gravitational force
F_grav_loc = weight*N_loc'*density*g*h;

%..Compute internal force
F_int_loc = weight*B_loc'*E*B_loc*h;

end