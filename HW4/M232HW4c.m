clc; close all; clear all;
%
% bvp_2.m 
% second order finite difference method for the bvp
%   u''(x) = f(x),   u'(ax)=sigma,   u(bx)=beta
% Using 3-pt differences on an arbitrary nonuniform grid.
% Should be 2nd order accurate if grid points vary smoothly, but may
% degenerate to "first order" on random or nonsmooth grids. 
%
% Different BCs can be specified by changing the first and/or last rows of 
% A and F.
%
% From  http://www.amath.washington.edu/~rjl/fdmbook/  (2007)

ax = 0;
bx = 1;
sigma = 3;   % Dirichlet boundary condition at ax
beta = -5;     %Neumann boundary condtion at bx

f = @(x) exp(x);  % right hand side function
utrue = @(x) exp(x) - (5+exp(1))*x + 2;%exp(x) + (sigma-exp(ax))*(x - bx) + beta - exp(bx);  % true soln

% true solution on fine grid for plotting:
xfine = linspace(ax,bx,101);
ufine = utrue(xfine);

% Solve the problem for ntest different grid sizes to test convergence:
m1vals = [10 20 40 80];
ntest = length(m1vals);
hvals = zeros(ntest,1);  % to hold h values
hvals2 = zeros(ntest,1);
E = zeros(ntest,1);      % to hold errors

for jtest=1:ntest
  m1 = m1vals(jtest);
  m2 = m1 + 1;
  m = m1 - 1;           % number of interior grid points
  hvals(jtest) = (bx-ax)/m1;  % average grid spacing, for convergence tests
  
  % Finer Grids (finding V)
  n1 = 2*m1vals(jtest);
  n2 = n1 + 1;                
  n = n1 - 1;
  hvals2(jtest) = (bx-ax)/n1;

  % set grid points:  
  gridchoice = 'uniform';          % see xgrid.m for other choices
  x = xgrid(ax,bx,m,gridchoice);  
  x2 = xgrid(ax,bx,n,gridchoice); %% FIner grid

  % set up matrix A (using sparse matrix storage):
  A = spalloc(m2,m2,3*m2);   % initialize to zero matrix
  A2 = spalloc(n2,n2,3*n2);

  % first row for Dirichlet BC at ax:
  A(1,1:3) = fdcoeffF(0, x(1), x(1:3)); 
  A2(1,1:3) = fdcoeffF(0, x(1), x(1:3)); 
  
  % interior rows:
  for i=2:m1
      A(i,i-1:i+1) = fdcoeffF(2, x(i), x((i-1):(i+1)));
  end
  
  for i=2:n1
      A2(i,i-1:i+1) = fdcoeffF(2, x2(i), x2((i-1):(i+1)));
  end

  % last row for Neumann BC at bx:
  A(m2,m-2:m2) = fdcoeffF(1,x(m2),x(m-2:m2)); 
  A2(n2,n-2:n2) = fdcoeffF(1,x2(n2),x2(n-2:n2)); 
  
  % Right hand side:
  F = f(x);
  F(1) = sigma;
  F(m2) = beta;
  
  F2 = f(x2);
  F2(1) = sigma;  
  F2(n2) = beta;
  
  % solve linear system:
  U = A\F;
  V = A2\F2;
  U(2:end-1) = (4*V(3:2:end-2)-U(2:end-1))/3;

  % compute error at grid points:
  uhat = utrue(x);
  err = U - uhat;
  err(end) = 0;
  E(jtest) = max(abs(err));  
  disp(' ')
  fprintf('Error with %i points is %9.5e\n',m2,E(jtest))

%   clf
  plot(x,U,'o')  % plot computed solution
  title(sprintf('Computed solution with %i grid points',m2));
  hold on
  plot(xfine,ufine)  % plot true solution
  hold off
  
  %pause to see this plot:  
  drawnow
  input('Hit <return> to continue ');
  
  end

error_table(hvals, E);   % print tables of errors and ratios
error_loglog(hvals, E);  % produce log-log plot of errors and least squares fit

