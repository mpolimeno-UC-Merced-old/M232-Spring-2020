t0 = 0;
T = 20;
alpha = 0.7;
beta = 0.7; 

t = linspace(t0,T,31);
thguess = @(x) 0.7+sin(x/2);  %initial guess %x is time here
m = length(t);
h = 1/20;
e = ones(m,1);
e = e*(h^-2);

thk = thguess(t)';
thkvec = [];
dkvec = [];
G = zeros(m,1);
Gvec = [];

for kk=1:10
    thkvec = [thkvec thk];
    G(1) = (1/h^2)*(alpha-2*thk(1)+thk(2)) + sin(thk(1));
    for pp = 2:m-1
        G(pp) = (1/h^2)*(thk(pp-1)-2*thk(pp)+thk(pp+1)) + sin(thk(pp));
    end
    G(m) = (1/h^2)*(thk(m-1)-2*thk(m)+beta) + sin(thk(m));
    Gvec = [Gvec G]; %to check if I am updating my G
    J = spdiags([e -2*e e],-1:1,m,m); %you need to reset your J at every iteration
    J = J + eye(length((J))).*cos(thk);
    delta = -J\G; %Newton steps
    thkp1 = thk + delta;
    dk = norm(delta,inf);
    dkvec = [dkvec dk];
    thk = thkp1; %update
    figure(1)
    plot(t,thk)
    hold on
    plot(t,thguess(t),'k--')
    title(sprintf('Computed solution after %i Newton Iterations',kk));
    xlim([t0 T])
    xlabel('$t$','interpreter','latex','fontsize',15)
    ylabel('$\theta^{[k]}$','interpreter','latex','fontsize',15)
    grid on
end

kvec = 1:kk;
figure(2)
loglog(kvec,dkvec,'b*-')
title(sprintf('Convergence for %i Newton Iterations',kk));
xlabel('$k$','interpreter','latex','fontsize',15)
ylabel('$\mid\mid{\delta^{[k]}}\mid\mid_{\infty}$','interpreter','latex','fontsize',15)
grid on