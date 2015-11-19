function [M_glob] = assemble_matrix(n_e, T, M_loc)

%..Assemble mass matrix
M_glob = T(:,:,1)'*M_loc*T(:,:,1);
for j=2:n_e
    Temp_m = T(:,:,j)'*M_loc*T(:,:,j);
    M_glob = M_glob+Temp_m;
    clear Temp_m
end

end