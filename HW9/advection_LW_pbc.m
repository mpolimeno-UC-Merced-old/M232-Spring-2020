function [h,k,error] = advection_LW_pbc(m,mlength)
%
% Solve u_t + au_x = 0  on [ax,bx] with periodic boundary conditions,
% using the Lax-Wendroff method with m interior points.
%
% Returns k, h, and the max-norm of the error.
% This routine can be embedded in a loop on m to test the accuracy,
% perhaps with calls to error_table and/or error_loglog.
%
% From  http://www.amath.washington.edu/~rjl/fdmbook/  (2007)

global a
for ii=1:length(mlength)

a = 2;           % advection velocity

clf              % clear graphics

ax = 0;
bx = 1;
tfinal = 0.1;                % final time

h = (bx-ax)/(m+1);         % h = delta x
k = 0.6*h ;                 % time step
nu = a*k/h;                % Courant number
x = linspace(ax,bx,m+2)';  % note x(1)=0 and x(m+2)=1
                           % With periodic BC's there are m+1 unknowns u(2:m+2)
I = 2:(m+2);   % indices of unknowns

nsteps = round(tfinal / k);    % number of time steps
%nplot = 20;       % plot solution every nplot time steps
                  % (set nplot=2 to plot every 2 time steps, etc.)
nplot = nsteps;  % only plot at final time

% if abs(k*nsteps - tfinal) > 1e-5
%    % The last step won't go exactly to tfinal.
%    disp(' ')
%    fprintf('WARNING *** k does not divide tfinal, k = %9.5e\n',k)
%    disp(' ')
% end

% initial conditions:
tn = 0;
u0 = eta(x);
u = u0;

% periodic boundary conditions:
u(1) = u(m+2);   % copy value from rightmost unknown to ghost cell on left
u(m+3) = u(2);   % copy value from leftmost unknown to ghost cell on right

% initial data on fine grid for plotting:
xfine = linspace(ax,bx,1001);
ufine = utrue(xfine,0);

% plot initial data:
figure(ii)
plot(x,u0,'b.-', xfine,ufine,'r')
axis([0 1 -.2 1.2])
legend('computed','true')
title(sprintf('Initial data at time = 0 with %5i grid points',m+1)) %m+1 beacuse at x(m+2)=x(1)
xlabel('$x$','interpreter','latex','fontsize',15)
ylabel('$u(x,t)$','interpreter','latex','fontsize',15)
%input('Hit <return> to continue  ');

% main time-stepping loop:

for n = 1:nsteps
     tnp = tn + k;   % = t_{n+1}

     % Lax-Wendroff:
     u(I) = u(I) - 0.5*nu*(u(I+1) - u(I-1)) + ...
                   0.5*nu^2 * (u(I-1) - 2*u(I) + u(I+1));

     % periodic boundary conditions:
     u(1) = u(m+2);   % copy value from rightmost unknown to ghost cell on left
     u(m+3) = u(2);   % copy value from leftmost unknown to ghost cell on right

     % plot results at desired times:
     %if mod(n,nsteps)==0 || n==nsteps
        uint = u(1:m+2);  % points on the interval (drop ghost cell on right)
        ufine = utrue(xfine,tnp);
        if n<nsteps
        figure(ii+1)
        plot(x,uint,'b.-', xfine,ufine,'r')
        axis([0 1 -.2 1.2])
        title(sprintf('t = %9.5e  after %4i time steps with %5i grid points',...
                       tnp,n,m+1))
        xlabel('$x$','interpreter','latex','fontsize',15)
        ylabel('$u(x,t)$','interpreter','latex','fontsize',15)
        error = max(abs(uint-utrue(x,tnp)));
        fprintf('at time t = %9.5e  max error =  %9.5e\n',tnp,error)
        end
%         if n<nsteps
%            input('Hit <return> to continue  '); 
%         end
     %end
     tn = tnp;   % for next time step
end
end
end
%--------------------------------------------------------

function utrue = utrue(x,t) %this gives the IC at t = dt
% true solution for comparison
global a

% For periodic BC's, we need the periodic extension of eta(x).
% Map x-a*t back to unit interval:

xat = rem(x - a*t, 1);
ineg = find(xat<0);
xat(ineg) = xat(ineg) + 1;
utrue = eta(xat);
return
end

%--------------------------------------------------------

function eta = eta(x) %this gives the IC at time 0
% initial data

beta = 600;
eta = exp(-beta*(x - 0.5).^2);
return
end

