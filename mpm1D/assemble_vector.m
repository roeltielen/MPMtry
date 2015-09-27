function [M] = assemble_mpm (n_n, n_e, elements_index, A)

%% Boolean matrices required for the assemblage procedure
T = zeros(2, n_n, n_e);
for j=1:n_e
   T(1,elements_index(j,1),j) = 1;
   T(2,elements_index(j,2),j) = 1;
end

%% Check whether A is a matrix or a vector
if ndims(A) == 2
    %assemble vector
    M = T(:,:,1)'*A(:,1);
    for i=2:n_e
        Temp = T(:,:,i)'*A(:,i);
        M = M + Temp;
        clear Temp_p
    end
 else    
    %assemble matrix
    M = T(:,:,1)'*A(:,:,1)*T(:,:,1);
    for j = 2:n_e
        Temp = T(:,:,j)'*A(:,:,j)*T(:,:,j);
        M = M+Temp;
        clear Temp
    end
    
end


