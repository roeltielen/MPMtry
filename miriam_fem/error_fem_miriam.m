function [] = error_fem_miriam()

[num_sol, exact_sol, M, time_step, dt] = fem_miriam();
 num_el = size(num_sol,1);
 x = 0:1/(num_el-1):1;
 plot(x,num_sol)
 hold on
 plot(x,exact_sol, '--r')

err    = num_sol - exact_sol;
errnrm = 0;
for i = 1:length(num_sol)
    errnrm = errnrm + M(i)*err(i)^2;
end
errnrm = sqrt(errnrm)

end