
% poisson2.m  -- solve the Poisson problem u_{xx} + u_{yy} = f(x,y)
% on [a,b] x [a,b].  
% 
% The 5-point Laplacian is used at interior grid points.
% This system of equations is then solved using backslash.
% 
% From  http://www.amath.washington.edu/~rjl/fdmbook/chapter3  (2007)
clc; clear all; close all;

ax = 0; 
bx = 1; 

ay = 0;
by = 2;

m1valx = [10 20 40 80];
m1valy = [21 41 81 161];
ntest = length(m1valx);
hvalx = zeros(ntest,1);
hvaly = zeros(ntest,1);
E = zeros(ntest,1);

for jtest = 1:ntest
    
    m = m1valx(jtest);
    my = m1valy(jtest);

    hx = (bx-ax)/(m+1);
    hy = (by-ay)/(my+1);
    x = linspace(ax,bx,m+2);   % grid points x including boundaries
    y = linspace(ay,by,my+2);   % grid points y including boundaries
    hvalx(jtest) = hx;
    hvaly(jtest) = hy;


    [X,Y] = meshgrid(x,y);      % 2d arrays of x,y values
    X = X';                     % transpose so that X(i,j),Y(i,j) are
    Y = Y';                     % coordinates of (i,j) point

    Iint = 2:m+1;              % indices of interior points in x
    Jint = 2:my+1;              % indices of interior points in y
    Xint = X(Iint,Jint);       % interior points
    Yint = Y(Iint,Jint);

    f = @(x,y) 1.25*exp(x+y/2);         % f(x,y) function

    rhs = f(Xint,Yint);        % evaluate f at interior points for right hand side
                               % rhs is modified below for boundary conditions.

    utrue = exp(X+Y/2);        % true solution for test problem

    % set boundary conditions around edges of usoln array:

    usoln = utrue;              % use true solution for this test problem
                                % This sets full array, but only boundary values
                                % are used below.  For a problem where utrue
                                % is not known, would have to set each edge of
                                % usoln to the desired Dirichlet boundary values.


    % adjust the rhs to include boundary terms:
    rhs(:,1) = rhs(:,1) - usoln(Iint,1)/hy^2;
    rhs(:,my) = rhs(:,my) - usoln(Iint,my+2)/hy^2;
    rhs(1,:) = rhs(1,:) - usoln(1,Jint)/hx^2;
    rhs(m,:) = rhs(m,:) - usoln(m+2,Jint)/hx^2;


    % convert the 2d grid function rhs into a column vector for rhs of system:
    F = reshape(rhs,m*my,1);

    % form matrix A:
    I = speye(my);
    I2 = speye(m)/hy^2;
    e = ones(m,1);
    e2 = ones(my,1);
    T = spdiags([e/hx^2 (-2/hx^2-2/hy^2)*e e/hx^2],[-1 0 1],m,m);
    S = spdiags([e2 e2],[-1 1],my,my);
    A = (kron(I,T) + kron(S,I2));


    % Solve the linear system:
    uvec = A\F;  

    % reshape vector solution uvec as a grid function and 
    % insert this interior solution into usoln for plotting purposes:
    % (recall boundary conditions in usoln are already set) 

    usoln(Iint,Jint) = reshape(uvec,m,my);

    % assuming true solution is known and stored in utrue:
    err = max(max(abs(usoln-utrue)));   
    fprintf('Error relative to true solution of PDE = %10.3e \n',err)

    % plot results:

    %err = max(abs(uvec - usoln));
    
    E(jtest) = err;
end

clf
hold on

figure(1);
contour(X,Y,utrue,30,'k')
hold on
contour(X,Y,usoln,30,'--r')
axis([ax bx ay by])
%daspect([1 1 1])
xlabel('x','interpreter','latex','fontsize',15)
ylabel('y','interpreter','latex','fontsize',15)
title('Contour Plot','interpreter','latex','fontsize',15)
%print('Contour1c','-depsc','-tiff')
legend('Analytical Solution','Numerical Solution','interpreter','latex')
hold off

figure(2)
surf(X,Y,usoln)
xlabel('x','interpreter','latex','fontsize',15)
ylabel('y','interpreter','latex','fontsize',15)
zlabel('u','interpreter','latex','fontsize',15)
colorbar
title('Numerical solution','interpreter','latex','fontsize',15)
%print('surfacenum1c','-depsc','-tiff')

figure(3)
surf(X,Y,utrue)
xlabel('x','interpreter','latex','fontsize',15)
ylabel('y','interpreter','latex','fontsize',15)
zlabel('u','interpreter','latex','fontsize',15)
colorbar
title('Analytical solution','interpreter','latex','fontsize',15)
%print('surfacetrue1c','-depsc','-tiff')

figure(4);
error_table(hvalx, E);   % print tables of errors and ratios
error_loglog(hvalx, E);

