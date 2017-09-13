function theta = angle_2ves(x1,y1,x2,y2)
% ANGLE_2VES
%   Find relative angle formed by the vector P1-P2 from horizontal line.
%   P1 = (a1, b1), P2 = (a2, b2) are the centroids of two vesicles.
%

[a1, b1] = centroid(x1,y1);
[a2, b2] = centroid(x2,y2);
theta = atan2d(b2-b1, a2-a1);


function [xbar, ybar] = centroid(x,y)
% CENTROID
%   This function computes the centroid (xbar,ybar) of a 2D region whose
%   boundary is given by (x(t),y(t)), parametrized by 
%   t = linspace(0,2*pi,n+1).
%

if nargin == 0
    n = 200;
    t = linspace(0,2*pi,n+1).';
    r = 8-sin(t)+2*sin(3*t)+2*sin(5*t)-sin(7*t)+3*cos(2*t)-2*cos(4*t);
    x = cos(t).*r;
    y = sin(t).*r;
    [xbar, ybar] = centroid(x,y);
    plot(x,y,'r',xbar,ybar,'o','linesmoothing','on')
    axis equal
    axis([-5,5,-5,5]*3)
    legend('butterfly', 'centroid')
end

n = length(x)-1;
xp = Dtp(x);
yp = Dtp(y);
A = pi/n*sum(x(2:end).*yp(2:end)-y(2:end).*xp(2:end));
xbar = pi/n/A*sum(x(2:end).^2.*yp(2:end));
ybar = -pi/n/A*sum(y(2:end).^2.*xp(2:end));

function dxdt=Dtp(x,d)

persistent K
%
% This program computes the derivative of a function of the form 
% f(t)=a*t+p(t), where a is a constant and p(t) is a periodic function on 
% the interval [0,2*pi]. Both endpoints must be included.
%
% Example:
%   t=linspace(0,2*pi,33)';
%   y=sin(t);
%   dy=Dt(y);
%

[m,p]=size(x);
m=m-1;
if isempty(K)||size(K,1)~=m||size(K,2)~=p
    if (-1)^m==1
        K=1i*[(0:m/2-1)';0;(-m/2+1:-1)']*ones(1,p);
    else
        K=1i*[(1:(m-1)/2)';(-(m-1)/2:0)']*ones(1,p);
    end
end

if nargin==1
    dxdt=ifft(K.*fft(x(1:m,:)));
else
    dxdt=ifft(K.^d.*fft(x(1:m,:)));
end
dxdt=real([dxdt;dxdt(1,:)]);