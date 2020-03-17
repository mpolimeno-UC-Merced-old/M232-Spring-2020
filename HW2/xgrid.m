function x = xgrid(ax,bx,m,gridchoice)
%
% Specify grid points in space for solving a 2-point boundary value problem
% or time-dependent PDE in one space dimension.
%
% Grid has m interior points on the interval [ax, bx].  
% gridchoice specifies the type of grid (see below).
% m+2 grid points (including boundaries) are returned in x.
%
% From  http://www.amath.washington.edu/~rjl/fdmbook/  (2007)

z = linspace(0, 1, m+2)';   % uniform grid in [0,1]

switch gridchoice

  case 'uniform'
     x = ax + (bx-ax)*z;

  case 'rtlayer'
     % Clustered near right boundary
     x =  ax + (bx-ax) * (1 - (1-z).^2);
     
  case 'random'
     x = ax + (bx-ax)*sort(rand(m+2,1));  
     x(1) = ax;  
     x(m+2) = bx;

  case 'chebyshev'
     % Chebyshev extreme points
     x = ax + (bx-ax) * 0.5*(1 + cos(pi*(1-z)));

  case 'legendre'
     % zeros of Legendre polynomial plus endpoints
     Toff = .5./sqrt(1-(2*(1:m-1)).^(-2));
     T = diag(Toff,1) + diag(Toff,-1);
     xi = sort(eig(T));
     x = ax + (bx-ax) * 0.5*(1 + [-1; xi; 1]);

  end
