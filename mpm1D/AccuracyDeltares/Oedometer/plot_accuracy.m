num_el = [4 8 16 32 64 128];
%accuracy FEM
error_fem = [0.138453265858977 0.039782741101408 0.010653172885301...
    0.003407251528314 0.002278295281541 0.002298645122159];
rel_error_fem = [0.010244555639088 0.002964810441508 7.953563062158936e-004...
    2.544969284958956e-004 1.701912694394808e-004  1.717162607007559e-004];
accuracy_fem = log2(error_fem(1:end-1)./error_fem(2:end))
accuracy_fem_rel = log2(rel_error_fem(1:end-1)./rel_error_fem(2:end))
%Richardson for FEM: node on top
N_4 = 23,77129786;
N_8 = 2,37712975E+01;
N_16 = 23,77129750;
N_32 = 23,77129750;
N_64 = 23,77129750;
N_128 = 23,77129750;

%accuracy FEM with load applied gradually
error_fem_gl = [0.730763226543127 0.718632972163952 0.714263190726859...
    0.713141503239686 0.712864241009948 0.712794255869300];
rel_error_fem_gl = [0.051657858603143 0.051123340469865 0.050893544533600...
    0.050833910910507 0.050819220178854 0.050815499726873
];
accuracy_fem_gl = log2(error_fem_gl(1:end-1)./error_fem_gl(2:end))
accuracy_fem_rel_gl = log2(rel_error_fem_gl(1:end-1)./...
    rel_error_fem_gl(2:end))

%all h (except 1/16) obviously provide the order of accuracy equal to 0

%accuracy MPM(4)
error_mpm4 = [0.235507443882444 0.100296213771872 0.035649697804837...
    0.016969106333397 0.012941060188519 0.011707256637032];
rel_error_mpm_4 = [0.017533259220214 0.007486051709825 0.002662509284523
0.001267544151858 9.666980339626251e-004 8.745459159557817e-004];
accuracy_mpm4 = log2(error_mpm4(1:end-1)./error_mpm4(2:end))
accuracy_mpm4_rel = log2(rel_error_mpm_4(1:end-1)./rel_error_mpm_4(2:end))




%accuracy matlab fem small deformations
%0.003019567353062 0.001870014408140 0.001704731044996 8.666809270091558e-004