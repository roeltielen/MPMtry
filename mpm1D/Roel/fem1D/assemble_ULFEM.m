function [M, F_grav, F_int_part] = assemble_ULFEM(n_e, T,  M_loc,...
 F_grav_loc, F_int_loc)

%..Assemble global mass matrix
%M = T(:,:,1)'*M_loc(:,:)*T(:,:,1);
M = T(:,:,1)'*M_loc(:,:,1)*T(:,:,1);
for j=2:n_e
   Temp_m = T(:,:,j)'*M_loc(:,:,j)*T(:,:,j);
   M = M+Temp_m;
   clear Temp_m
end

%..Assemble global gravitational force vector
F_grav = T(:,:,1)'*F_grav_loc(:,:,1);
for k=2:n_e
    Temp_g = T(:,:,k)'*F_grav_loc(:,:,k);
    F_grav = F_grav+Temp_g;
    clear Temp_g
end

%..Assemble global internal force vector 
F_int_part = T(:,:,1)'*F_int_loc(:,:,1)*T(:,:,1);
for l=2:n_e
   Temp_i = T(:,:,l)'*F_int_loc(:,:,j)*T(:,:,l);
   F_int_part = F_int_part+Temp_i;
   clear Temp_i
end



