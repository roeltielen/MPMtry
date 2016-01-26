num_el = [4 8 16 32 64 128];
%accuracy MPM(1)
% Deltares
error_1 = [0.002567737343201 8.194641928163814e-004 2.259255963095032e-004...
    7.140661483842501e-005 2.050956700532768e-005 3.456171713329815e-006];
loglog(num_el, error_1)
hold on
accuracy_1 = log2(error_1(1:end-1)./error_1(2:end))
% Matlab
error_matlab_1 =  [1.2918E-03 3.2595E-04 8.1795E-05 2.0632E-05...
    5.4969E-06 2.1502E-06];
accuracy_matlab_1 = log2(error_matlab_1(1:end-1)./error_matlab_1(2:end))

%accuracy MPM(4)
error_4 = [0.001456578325061 4.149949874748430e-004 1.062060321475216e-004...
    2.895357726837570e-005 9.583630043251487e-006 4.955854241281396e-006];
accuracy_4 = log2(error_4(1:end-1)./error_4(2:end))
loglog(num_el, error_4,'r')
hold on
error_matlab_4 = [1.0374e-03 2.6167e-04 6.5694e-05 1.6625e-05 4.5505e-06...
    1.9769e-06];
accuracy_matlab_4 = log2(error_matlab_4(1:end-1)./error_matlab_4(2:end))

%accuracy MPM(8)
error_8 = [0.001697435797358 0.034818575189317 1.291696359957010e-004...
        3.701819141961082e-005 1.147304338790864e-005 4.311498131850748e-006];
accuracy_8 = log2(error_8(1:end-1)./error_8(2:end))
loglog(num_el, error_8,'k')