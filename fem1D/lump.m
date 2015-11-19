function[M_lumped] = lump(M)
    M_diag = sum(M,2);
    M_lumped = diag(M_diag);    
end