function [errnrm] = compute_error_norm(num_sol, exact_sol, M_lump)
% Computes the norm of the discretization error in 2-norm

err    = num_sol - exact_sol;
errnrm = 0;
for i = 1:length(num_sol)
    errnrm = errnrm + M_lump(i,i)*err(i)^2;
end
errnrm = sqrt(errnrm);