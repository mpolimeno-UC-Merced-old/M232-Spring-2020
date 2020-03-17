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
bx = pi;
alpha = linspace(-2,2,11);   % Dirichlet boundary condition at ax
beta = -alpha;     % Dirichlet boundary condtion at bx

f = @(x) 0;  % right hand side function
%utrue = @(x) sin(x)+alpha*cos(x);  % true soln

% true solution on fine grid for plotting:
xfine = linspace(ax,bx,101);
%ufine = utrue(xfine);

% Solve the problem for ntest different grid sizes to test convergence:
m1vals = [10 20 40 80];
ntest = length(m1vals);
hvals = zeros(ntest,1);  % to hold h values
E = zeros(ntest,1);      % to hold errors

for jtest=1:ntest
  m1 = m1vals(jtest);
  m2 = m1 + 1;
  m = m1 - 1;                 % number of interior grid points
  hvals(jtest) = (bx-ax)/m1;  % average grid spacing, for convergence tests

  % set grid points:  
  gridchoice = 'uniform';          % see xgrid.m for other choices
  x = xgrid(ax,bx,m,gridchoice);   %this is what I am feeding to my f(x)

  % set up matrix A (using sparse matrix storage):
  A = spalloc(m2,m2,3*m2);   % initialize to zero matrix (3*m2 gives the number of nonzero elements. But why 3*m2??)
                            %probably b/c it uses a 3-point difference
  % first row for Dirichlet BC at ax:
  A(1,1:3) = fdcoeffF(0, x(1), x(1:3)); %k=0 is the zeroth order derivative

  % interior rows:
  for i=2:m1
     A(i,i-1:i+1) = fdcoeffF(2, x(i), x((i-1):(i+1)));
  end
  
  for ii = 2:m1
      A(ii,ii) = A(ii,ii)+1;
  end

  % last row for Dirichlet BC at bx:
  A(m2,m:m2) = fdcoeffF(0,x(m2),x(m:m2)); %k=0 is the zeroth order derivative


  % compute error at grid points:
  for jj=1:length(alpha)
      utrue = @(x) sin(x)+alpha(jj)*cos(x);  % true soln
      ufine = utrue(xfine);
      plot(xfine,ufine)  % plot true solution
      hold on
      grid on
      xlim([xfine(1) xfine(end)])
  end
  
  Legend=cell(2,1);
  for kk=1:length(alpha)
      Legend{kk} = strcat('\alpha = ',num2str(alpha(kk)));
  end
  legend(Legend,'location','best')
  xlabel('$x$','interpreter','latex','fontsize',15)
  ylabel('$u(x)$','interpreter','latex','fontsize',15)
  title('True Solution for different $\alpha$ values','interpreter','latex','fontsize',15)
  
end

