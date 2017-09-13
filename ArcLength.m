function A=VWM_ArcLength(X,Y,l)
% This function computes the arc length of a vesicle with parameterization
% X(t) and Y(t). t must be define on an interval [0,l), where l is the
% width of the interval.
%
% Example:
%   t=linspace(0,2*pi,65)';
%   x=cos(t);
%   y=sin(t);
%   VWM_ArcLength(x,y,2*pi)
%

m=size(X,1);
sa=sqrt(Dtp(X).^2 + Dtp(Y).^2);
A=l/(m-1)*sum(sa(1:m-1,:))';