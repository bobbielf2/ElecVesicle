function z=INTt(x)

% This function computes the integral of p(t) from t=0 to t=x, where p(t) 
% is a periodic function on [0,2*pi].
%
% Example:
%   t=linspace(0,2*pi,33)';
%   y=sin(t);
%   Iy=INTt(y);
%

[m,p]=size(x);
t=linspace(0,2*pi,m)';
m=m-1;
K=1./((-floor(m/2):ceil(m/2)-1)'*1i*ones(1,p));
if (-1)^m==1
K(1)=0;
end
K(ceil(m/2)+1)=0;
b=sum(x(1:m,:),1)/m;
z=ifft(ifftshift(K.*fftshift(fft(x(1:m,:)))));
z=real([z;z(1,:)])+t*b;
z=z-z(1);
end