%function [h,k,error] = heat_CN(m)
%
% heat_CN.m
%
% Solve u_t = kappa * u_{xx} on [ax,bx] with Dirichlet boundary conditions,
% using the Crank-Nicolson method with m interior points.
%
% Returns k, h, and the max-norm of the error.
% This routine can be embedded in a loop on m to test the accuracy,
% perhaps with calls to error_table and/or error_loglog.
%
% From  http://www.amath.washington.edu/~rjl/fdmbook/  (2007)

clf              % clear graphics
%hold on          % Put all plots on the same graph (comment out if desired)

ax = 0;
bx = 1;
kappa = .02;               % heat conduction coefficient:
tfinal = 1;                % final time
errvec = [];

mvals = [8 18 38 68 98];
hvals = [];
kvals = [];
for ii=1:length(mvals)
    m = mvals(ii);
    h = (bx-ax)/(m+1);         % h = delta x
    hvals = [hvals h];
    x = linspace(ax,bx,m+2)';  % note x(1)=0 and x(m+2)=1
                           % u(1)=g0 and u(m+2)=g1 are known from BC's
    k = 4*h; %just picking a k=O(h)                  % time step
    kvals = [kvals k];
    
    
    nsteps = round(tfinal / k);    % number of time steps
    %nplot = 1;      % plot solution every nplot time steps
                 % (set nplot=2 to plot every 2 time steps, etc.)
    nplot = nsteps;  % only plot at final time

%     if abs(k*nsteps - tfinal) > 1e-5
%         % The last step won't go exactly to tfinal.
%         disp(' ')
%         fprintf('WARNING *** k does not divide tfinal, k = %9.5e\n',k)
%         disp(' ')
%     end


    % true solution for comparison:
    % For Gaussian initial conditions u(x,0) = exp(-beta * (x-0.4)^2)
    beta = 150;
    utrue = @(x,t) exp(-(x-0.4).^2 / (4*kappa*t + 1/beta)) / sqrt(4*beta*kappa*t+1);

    % initial conditions:
    u0 = utrue(x,0);


    % Each time step we solve MOL system U' = AU + g using the Trapezoidal method

    % set up matrices:
    r = kappa* k/(4*h^2);
    e = ones(m,1);
    A = spdiags([e -2*e e], -1:1, m, m);
    A1 = eye(m) - r * A;
    A2 = eye(m) + r * A;
    
    %matrix stuff to go from nphalf to np1
    r2 = kappa*k/(3*h^2);
    A3 = eye(m) - r2*A;
    
    % initial data on fine grid for plotting:
    xfine = linspace(ax,bx,1001);
    ufine = utrue(xfine,0);

    % initialize u and plot:
    tn = 0;
    u = u0;
    
    figure(ii)
    plot(x,u,'b.-', xfine,ufine,'r')
    legend('TR-BDF2','Analytical','interpreter','latex')
    title(sprintf('Initial data at $t = 0$ for %5i grid points',m+2),'interpreter','latex','fontsize',15)
    xlabel('$x$','interpreter','latex','fontsize',15)
    ylabel('$u(x,t)$','interpreter','latex','fontsize',15)

    %input('Hit <return> to continue  ');


    % main time-stepping loop:

    for n = 1:nsteps
        tnp = tn + k;   % = t_{n+1} %total number of time steps
        tnstar = tn+(k/2); % = t_{n+1/2}
        % boundary values u(0,t) and u(1,t) at times tn and tnp:
        g0n = u(1);
        g1n = u(m+2);
        g0np = utrue(ax,tnp);
        g0star = utrue(ax,tnstar);
        g1star = utrue(bx,tnstar);
        g1np = utrue(bx,tnp);

        % compute right hand side for linear system:
        uint = u(2:(m+1));   % interior points (unknowns)
        rhsstar = A2*uint;
        % fix-up right hand side using BC's from n to nphalf
        rhsstar(1) = rhsstar(1) + r*(g0n + g0star);
        rhsstar(m) = rhsstar(m) + r*(g1n + g1star);

        % solve linear system:
        uintstar = A1\rhsstar;
        
        rhsstep = (1/3)*(4*uintstar - uint); %this come from setting up tr-bdf2 and plugging in unphalf
        rhsstep(1) = rhsstep(1) + r2*(g0np);
        rhsstep(m) = rhsstep(m) + r2*(g1np);
        
        uint = A3\rhsstep; 
        % augment with boundary values:
        u = [g0np; uint; g1np];
        % plot results at desired times:
        if mod(n,nplot)==0 || n==nsteps
            ufine = utrue(xfine,tnp);
            figure(ii+5)
            plot(x,u,'b.-', xfine,ufine,'r')
            legend('Crank-Nicolson','Analytical','interpreter','latex')
            title(sprintf('t = %9.5e  after %4i time steps with %5i grid points',...
                       tnp,n,m+2),'interpreter','latex')
            xlabel('$x$','interpreter','latex','fontsize',15)
            ylabel('$u(x,t)$','interpreter','latex','fontsize',15)
            error = max(abs(u-utrue(x,tnp)));
            errvec = [errvec error];
            fprintf('at time $t$ = %9.5e  max error =  %9.5e\n',tnp,error)
%             if n<nsteps, input('Hit <return> to continue  '); 
%             end
        end

     tn = tnp;   % for next time step
    end
end

error_table(hvals, kvals, errvec);   % print tables of errors and ratios

figure(ii+6)
error_loglog(hvals, errvec);  % produce log-log plot of errors and least squares fit
