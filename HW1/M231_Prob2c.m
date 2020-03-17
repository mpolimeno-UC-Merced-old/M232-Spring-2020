%Problem 2c - HW1 for M232, January 2020

uptrue = -4*sin(2);
hvals = logspace(-1,-4,13);
Eupp = [];  E3u = [];
u = @(x) sin(2*x);
c = [-1/12;4/3;-5/2;4/3;-1/12]; %coefficients obtained in Problem 2b
% table headings:
disp(' ')
disp('       h          upp             D3u')


for i=1:length(hvals)
   h = hvals(i);
   % approximations to u''(1):
   %this approximation is as Problem 2b
   Dm2u = (h^-2)*c(1)*u(1-2*h);
   Dmu = (h^-2)*c(2)*u(1-h);
   D0u = (h^-2)*c(3)*u(1);
   Dpu = (h^-2)*c(4)*u(1+h);
   Dp2u = (h^-2)*c(5)*u(1+2*h);
   upp = Dm2u+Dmu+D0u+Dpu+Dp2u;
   
   %this is to compute coefficients with fdcoeffF and then compare with
   %what we found in Problem 2b
   xpts = 1 + h*(-2:2)';
   D3u = fdcoeffF(2,1,xpts)*sin(2*xpts) ;
   
   % errors:
   Eupp(i) = upp - uptrue;
   E3u(i) = D3u - uptrue;

   % print line of table:
   fprintf('%13.4e   %13.4e  %13.4e\n',...
                 h,Eupp(i),E3u(i))
end

y = hvals.^(4); %get line of slope 4 to check order of the method
% plot errors:
clf
loglog(hvals,y,'k--')
hold on
loglog(hvals,abs(Eupp),'ro-')
hold on
loglog(hvals,abs(E3u),'bo-')
grid on
legend('${\rm Slope}~4$','${\rm Vandemonde}$','{\rm FD-Coeffs}','interpreter','latex','location','best')
title('Error for $u^{\prime\prime}(x)$ for $u(x)=\sin(2x)$','interpreter','latex','fontsize',15)
xlabel('$h$','interpreter','latex','fontsize',15)
ylabel('$\mid{u(\bar{x})-u_{\rm true}(x)}\mid$','interpreter','latex','fontsize',15)

