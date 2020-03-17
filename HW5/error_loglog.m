function error_loglog(h,E)
%
% Produce log-log plot of E vs. h.
% Estimate order of accuracy by doing a linear least squares fit.
%
% From  http://www.amath.washington.edu/~rjl/fdmbook/  (2007)

h = h(:);            % make sure it's a column vector
E = E(:);            % make sure it's a column vector
ntest = length(h);
clf
loglog(h,E,'bo-')
axis([.5*min(h) 1.5*max(h)  .5*min(E) 1.5*max(E)])
title('log-log plot of errors vs. h','interpreter','latex','fontsize',20)
xlabel('$h$','interpreter','latex','fontsize',15)
ylabel('$\mid{{\bf U}-\hat{\bf U}}\mid$','interpreter','latex','fontsize',15)

% Estimate order of accuracy from least squares fit:
Ap = ones(ntest,2);
Ap(:,2) = log(h);
bp = log(E);
Kp = Ap\bp;
K = Kp(1);
p = Kp(2);
disp(' ')
fprintf('Least squares fit gives E(h) = %g * h^%g\n',exp(K),p)
disp(' ')

% add graph of this line to loglog plot:
hold on
err1 = exp(K)*h.^p;
loglog(h,err1,'r')
legend('errors', 'least squares fit','Location','SouthEast')
hold off
grid on
