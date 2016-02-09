%..to run this code on a single mesh, do n = <value>; poisson_2d 
%..to run thuvecis code for convergence analysis, do e.g. 
%..hvec = []; errnrmvec = []; 
%..for p=2:8 n = 2^p; poisson_1d; hvec = [hvec,1/n]; errnrmvec = [errnrmvec;errnrm]; end
%..loglog(hvec, errnrmvec,'*-')

%..Set number of elements and number of nodes (inclusing the boundary nodes)
%n = 3;
n1 = n-1; 
nnodes = (n+1)^2; 

%..Generate mesh size and vector of grid points.. 
h     = 1/n; 
h2    = h*h; 
xvec  = [0:h:1];
[X,Y] = meshgrid(xvec,xvec);  

%..Find global indices of interior and boundary nodes..
indic  = ones(n+1,n+1); indic(2:end-1,2:end-1) = 0
indvec = reshape(indic,nnodes,1); 
bnd    = find(indvec==1); interior = find(indvec==0);

%..Set source function and exact solution..
syms xs ys
exacts = xs^2*cos(pi*xs+pi*ys);
srcs   = -diff(diff(exacts,xs),xs)-diff(diff(exacts,ys),ys); 
src    = matlabFunction(srcs); 
%exact    = inline('1/6*y^3','x','y');  
exact  = matlabFunction(exacts);
%src  = inline('-y','x','y'); 

%..Construct 1D identity and the 1D system matrix
e = ones(n1,1); 
I1 = spdiags(e,0,n1,n1);
A1 = spdiags([-e 2*e -e],-1:1,n1,n1);
A1 = (1/h2)*A1; 

%..Construct the 2D system matrix by Kronecker products
A2 = kron(I1,A1) + kron(A1,I1);
A  = speye(nnodes,nnodes); 
A(interior,interior) = A2; 

%..Construct the exact solution and the rhs as 2-by-2 arrays 
ue = exact(X,Y);
f  = src(X,Y); 
%....Set correct values on boundary in the rhs-array
f(1,:) = ue(1,:); f(end,:) = ue(end,:); 
f(:,1) = ue(:,1); f(:,end) = ue(:,end); 
%....Set correct values in nodes connected to the boundary in the rhs-vector
f(2:end-1,2)     = f(2:end-1,2)     + (1/h2)*f(2:end-1,1); 
f(2:end-1,end-1) = f(2:end-1,end-1) + (1/h2)*f(2:end-1,end); 
f(2,2:end-1)     = f(2,2:end-1)     + (1/h2)*f(1,2:end-1); 
f(end-1,2:end-1) = f(end-1,2:end-1) + (1/h2)*f(end,2:end-1); 

%..Reshape the rhs-array to a rhs-vector 
fvec = reshape(f',nnodes,1);

%..Print the matrix and the rhs
full(A)
full(fvec)

% %..Solution
% %..Solve the discrete system of equations.. 
% uhvec = A \ fvec; 
% 
% %..Reshape the discrete solution to an array 
% uh = reshape(uhvec,n+1,n+1); uh = uh'; 
% 
% %..Plot the finite element solution superimposed on the exact solution.. 
% imagesc(ue-uh), colorbar
%  
% %..compute the norm of the discretization error in max-norm 
% err    = ue - uh; 
% %errnrm = norm(err,'inf'); 
% errvec = reshape(err,(n+1)^2,1);
% errnrm = norm(errvec, 'inf');
% fprintf('  The discretization errnrm = %e\n', errnrm)
