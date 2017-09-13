function [val, Ar, Le] = reduced_volumex(X)
% Given the radius at different time steps, compute the corresponding
% reduced volumes

[n, m] = size(X); n=n/2;

c0 = fftshift(1i*(-n/2:1:n/2-1)');
D1FT(c0, c0);

for j = 1:m
    X1 = X(1:n,j); X2 = X(1+n:2*n,j); 
    X1a = D1FT(X1); X2a = D1FT(X2);
    L = 2*pi/n*sum(sqrt(X1a.^2 + X2a.^2));
    A = pi/n*sum(X1.*X2a - X2.*X1a);
    Ar(j) = A;
    Le(j) = L;
    val(j) = A/(L^2/4/pi);
end

function val = D1FT(f, kIn)
% Computes m^{th}  derivative of a periodic function

persistent k  
if nargin>1
    k = kIn;
end

val = real(ifft(k.*fft(f)));
end

end
