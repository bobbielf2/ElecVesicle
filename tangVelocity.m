function val = tangVelocity(X,Vel)

n = length(X)/2;
c0 = fftshift(1i*(-n/2:1:n/2-1)'); D1 = @(f) real(ifft(c0.*fft(f)));
X1 = X(1:n); X2 = X(n+1:2*n);
sa = sqrt(D1(X1).^2 + D1(X2).^2);
tang = [D1(X1)./sa, D1(X2)./sa];

val = dot(tang,Vel,2);


% plot(X1, X2); hold on; axis equal
% quiver(X1,X2,tang(:,1),tang(:,2));