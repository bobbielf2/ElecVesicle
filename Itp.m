function xout=Itp(x,p)
%
% This function interpolates a function of the form f(t)=a*t+q(t) at the
% points p on the interval [0,2*pi]. q(t) must be a periodic function.
%
% Example:
%   t=linspace(0,2*pi,33)';
%   s=2*pi*rand(100,1);
%   It(sin(t),s)-sin(s)
%

m=size(x,1)-1;
X=fft(x(1:m,:));
K=[0:ceil(m/2)-1,-floor(m/2):-1];
xout=real(exp(1i*p*K)*X/m);