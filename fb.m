function X=fb(x,y,dsdt,kb)
%
% This function computes the bending force on a vesicle. 'x' and 'y' are
% vectors that contain the x and y-coordinates of the vesicle. 'dsdt' is a
% vector that contains the derivative of arc length with respect to the
% parameterization variable. 'dsdt' is computed using Dsdt.m. 'kb' is the
% bending modulus.
%

% We compute the 4th derivative of x and y with respect to arc length.
for k=1:4
DT=Dtp([x,y]);
x=DT(:,1)./dsdt;
y=DT(:,2)./dsdt;
end

% We multiply by the bending modulus.
X=-kb*[x,y];
end