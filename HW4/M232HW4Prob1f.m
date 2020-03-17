%
% -------------------------
% MATLAB: SpectralMethod1.m
% -------------------------
%
% This matlab code computes the solution to the Dirichlet, two-point
% boundary value problem of the form
% 
% u'' = exp(x) on [0,1], u(0) = 3, u(1) = -5
% 
% using a spectral collocation method on a uniform grid.
%

clear all;

% set the values of the Dirichlet and Neumann boundary conditions

alpha = 3.0;
beta  =  -5.0;

ax = 0;
bx = 1;

% set the number of interior grid points

mvals = 1:4:100;
for ii = 1:length(mvals)
 m = mvals(ii);
% compute the mesh width

h = (bx-ax) / ( m + 1.0 );

% compute the grid

%x = 0 : h : 3.0;

gridchoice = 'chebyshev';          % see xgrid.m for other choices
x = xgrid(ax,bx,m,gridchoice);

% compute the spectral collocation matrix

A = zeros( m + 2 );

% set the Dirichlet BC at x = 0

A(1,1) = 1.0;

% compute the interior rows

for j = 2 : m + 1
    
    A(j,:) = fdcoeffF( 2, x(j), x );
    
end;

% set the Neumann BC at x = 3

A(m+2,1:m+2) = fdcoeffF(1,x(m+2),x);

% compute the right-hand side

f      = zeros( m + 2, 1 );
f      = exp( x );
f(1)   = alpha;
f(m+2) = beta;

% solve the linear system of equations

u = A \ f;

% compute the exact solution

%C0 = alpha - 1.0;
%C1 = ( beta - alpha + 1.0 - exp( 3.0 ) ) / 3.0;

%uexact = C0 + C1 * x + exp( x );
uexact = exp(x) - (5+exp(1))*x + 2;

error_inf(ii) = max( abs( uexact - u ) );
error_one(ii) = h * sum( abs( uexact - u ) );
error_two(ii) = sqrt( h * sum( abs( uexact - u ).^2 ) );
%plot the solution

figure(ii)
plot( x, u, 'ko-', x, uexact, 'r--' );
xlabel('x','interpreter','latex','fontsize',15 );
ylabel('u(x)','interpreter','latex','fontsize',15 );
grid on
title( sprintf('Spectral Collocation Method on a Chebyshev Grid for %i points',m),'interpreter','latex');
legend( 'Spectral Collocation','Exact','interpreter','latex')
end

figure(ii+1)
errgrid = 1:length(mvals);
loglog(errgrid,error_inf)
hold on
loglog(errgrid,error_one)
hold on
loglog(errgrid,error_two)
grid on
xlabel('m','interpreter','latex','fontsize',15)
ylabel('Error','interpreter','latex','fontsize',15)
title('Error in different norms','interpreter','latex','fontsize',15)
legend('\mid{E}\mid_{\infty}','\mid{E}\mid_{1}','\mid{E}\mid_{2}')