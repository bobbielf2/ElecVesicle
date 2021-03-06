function Z=VWM_Area(x,y)
%
% This function computes the area of a vesicle. 'x' contains the
% x-coordinates and 'y' contains the vesicle's y-coordinates. The first
% points is need at the beginning and end of the vector.
%
% Example:
%   t=linspace(0,2*pi,65)';
%   x=cos(t);
%   y=sin(t);
%   VWM_Area(x,y)
%

int=x.*Dtp(y);
m=length(int);
Z=sum(int(1:m-1))*2*pi/(m-1);