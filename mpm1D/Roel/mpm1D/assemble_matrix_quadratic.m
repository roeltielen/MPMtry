function [M] = assemble_matrix_quadratic(n_n, n_e, elements_index, A)

%% Boolean matrices required for the assemblage procedure
T = zeros(3, n_n, n_e);
for j=1:n_e
    T(1,elements_index(j,1),j) = 1;
    T(2,elements_index(j,2),j) = 1;
    T(3,elements_index(j,3),j) = 1;
end


%assemble matrix
M = T(:,:,1)'*A(:,:,1)*T(:,:,1);
for j = 2:n_e
    Temp = T(:,:,j)'*A(:,:,j)*T(:,:,j);
    M = M+Temp;
    clear Temp
end

end


