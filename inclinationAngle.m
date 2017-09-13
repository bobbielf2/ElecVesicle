function val = inclinationAngle(X)
% bugs:
%   1. the angle returned by this function is between -pi/4 and 3*pi/4
%   2. when the ellipse inc angle is exactly pi/2 (90 degrees), it returns -pi/2

if nargin == 0, test_inclinationAngle; return; end

n = length(X)/2;
c0 = fftshift(1i*(-n/2:1:n/2-1)'); D1 = @(f) real(ifft(c0.*fft(f)));
X1 = X(1:n); X2 = X(n+1:2*n);
sa = sqrt(D1(X1).^2 + D1(X2).^2);
nor = [D1(X2)./sa, -D1(X1)./sa];

c = [mean(X1), mean(X2)];
r = [X1-c(1), X2-c(2)];
rho2 = r(:,1).^2 + r(:,2).^2;
rn = sum(r.*nor, 2);

intg = @(f) 1/4*2*pi/n*sum(f.*sa.*rn);
J(1,1) = intg(rho2 - r(:,1).^2);
J(2,2) = intg(rho2 - r(:,2).^2);
J(1,2) = intg(-r(:,1).*r(:,2)); 
J(2,1) = J(1,2);

[V,D] = eig(J);
[~, ind] = min(D);

v = V(:,ind(1));
val = -atan2(v(1), v(2));

function test_inclinationAngle

t = linspace(-pi,pi,65).'; t(end)=[];
x = 1.5*cos(t); y = 1*sin(t);

for th = linspace(0,2*pi,200)
  % th = 0;
  A = [cos(th),-sin(th);sin(th),cos(th)];
  xx = [x,y]*A(1,:).';
  yy = [x,y]*A(2,:).';
  ang = inclinationAngle([xx;yy])/pi*180;
  
  figure(1), clf
  plot(xx,yy,'r')
  axis equal
  axis([-2,2,-2,2])
  title(ang)
  drawnow
  pause(0.01)
end
% plot(X1, X2); hold on; rmax =sqrt(max(rho2)); plot([c(1), rmax*cos(val)], [c(2), rmax*sin(val)],'k-.');axis equal
% quiver(X1,X2,nor(:,1),nor(:,2));