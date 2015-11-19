function [F_glob] = assemble_vector(n_e, T, F_loc)

%..Assemble vector
F_glob = T(:,:,1)'*F_loc;
for j=2:n_e
    Temp_f = T(:,:,j)'*F_loc;
    F_glob = F_glob+Temp_f;
    clear Temp_f
end

end