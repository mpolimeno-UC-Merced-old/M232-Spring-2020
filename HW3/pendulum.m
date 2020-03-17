%
% bvp_2.m 
% second order finite difference method for the bvp
%   u''(t) = f(t),   u'(at)=sigma,   u(bt)=beta
% Using 3-pt differences on an arbitrary nonuniform grid.
% Should be 2nd order accurate if grid points vary smoothly, but may
% degenerate to "first order" on random or nonsmooth grids. 
%
% Different BCs can be specified by changing the first and/or last rows of 
% A and F.
%
% From  http://www.amath.washington.edu/~rjl/fdmbook/  (2007)

t0 = 0;
T = 2*pi;
alpha = 0.7;   % Neumann boundary condition at at
beta = 0.7;     % Dirichlet boundary condtion at bt

%g = @(theta) 0;  % right hand side function
thguess = @(x) 0.7*cos(x)+0.5*sin(x);  %initial guess %x is time here

% true solution on fine grid for plotting:
%tvec = linspace(t0,T,101);
%thfine = thguess(tvec);
%A = spdiags([e -2*e e],-1:1,n,n);
% Solve the problem for ntest different grid sizes to test convergence:
m1vals = [10 20 40 80];
ntest = length(m1vals);
hvals = zeros(ntest,1);  % to hold h values

for jtest=1:ntest
  m1 = m1vals(jtest);
  m2 = m1 + 1;
  m = m1 - 1;                 % number of interior grid points
  hvals(jtest) = (T-t0)/m1;  % average grid spacing, for convergence tests

  % set grid points:  
  gridchoice = 'uniform';          % see tgrid.m for other choices
  t = xgrid(t0,T,m,gridchoice); %this is my time variable
  
  %thk = thguess(t(1));
  G = zeros(m2,1);
  %J = spalloc(m2,m2,3*m2);
  for kk=1:8
    for pp = 2:m1
        G(1) = hvals(jtest).^(-2)*(alpha + thguess(t(1))-2*thguess(t(2)) + thguess(t(3))) + sin(thguess(t(2)));
        G(m2) = hvals(jtest).^(-2)*(thguess(t(end-3))-2*thguess(t(end-2)) + beta) + sin(thguess(t(end-1)));
        G(pp) = hvals(jtest).^(-2)*(thguess(t(pp-1))-2*thguess(t(pp)) + thguess(t(pp+1))) + sin(thguess(t(pp)));
        J = hvals(jtest).^(-2)*spdiags([1 -2+hvals(jtest).^(2)*cos(thguess(pp-1)) 1],1:1,m2,m2);
        delta = -J\G;
        thk = thguess(t(pp-1)) + delta;
        thguess(t(pp)) = thk;
    end
  end


%   % last row for Dirichlet BC at bt:
%   J(m2,m:m2) = fdcoeffF(0,t(m2),t(m:m2));
% %   G(m2,m:m2) = fdcoeffF(0,x(m2),x(m:m2));
%   
%   % Right hand side:
%   G = g(thguess(jtest));
%   G(1) = alpha;  
%   G(m2) = beta;
%   
%   % solve nonlinear system: 
%   thk = thguess(t(1)); %initial guess evaluated at first grid point
  %J2 = @(theta) J;
  %G2 = @(theta) G;
  
%  for jj=2:m1 %this is the part that keeps bugging me
%      J(jj,jj) = J(jj,jj)+cos(thguess(x(jj)));
%      J2 = @(theta) J;
%      G2 = @(theta) G;
%      for kk = 1:8 %newton iteration
%         delta = -J\G'; %this way I am not updating the jacobian (problem)
%         thkp1 = thk + delta;
%         thk = thkp1;
%      end
%   end
  
  
  % compute error at grid points:
%   uhat = thguess(t);
%   err = U - uhat;
%   E(jtest) = mat(abs(err));  
%   disp(' ')
%   fprintf('Error with %i points is %9.5e\n',m2,E(jtest))

%   clf
  plot(t,thk)  % plot computed solution
  title(sprintf('Computed solution with %i grid points',m2));
  hold on
%   plot(tvec,thfine)  % plot true solution
%   hold off
  
  % pause to see this plot:  
  drawnow
  input('Hit <return> to continue ');
  
  end

% error_table(hvals, E);   % print tables of errors and ratios
% error_loglog(hvals, E);  % produce log-log plot of errors and least squares fit


