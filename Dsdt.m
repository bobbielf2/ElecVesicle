function dsdt=Dsdt(x,y)
%
% This function computes the derivative of arc length with respect to a
% parameterization t on [0,2*pi]. The x and y-coordinates of the curve are
% given by 'x' and 'y'. Both functions must be described using functions fo
% the form x(t)=a*t+p(t), where a is a constant and p(t) is a periodic
% funciton on [0,2*pi]. Both end points must be included.
%
% Example:
%   t=linspace(0,2*pi,65)';
%   x=2*sin(t);
%   y=cos(t);
%   dsdt=Dsdt(x,y);
%

DT=Dtp([x,y]);
dsdt=sqrt(DT(:,1).^2+DT(:,2).^2);
end