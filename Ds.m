function [dxds,dyds]=Ds(x,y,dsdt)
%
% This function computes the derivative of x and y with respect to arc
% length. x(t) and y(t) must be functions of the form x(t)=a*t+p(t), where
% a is a constant and p(t) is a periodic function on [0,2*pi]. Both
% end points must be included.
%
% Example:
%   t=linspace(0,2*pi,65)';
%   x=2*sin(t);
%   y=cos(t);
%   dsdt=Dsdt(x,y);
%   [dxds,dyds]=Ds(x,y,dsdt);

DT=Dtp([x,y]);
dxds=real(DT(:,1)./dsdt);
dyds=real(DT(:,2)./dsdt);
end