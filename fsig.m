function X=fsig(x,y,sigma,dsdt)
%
% This function computes the bending force on a vesicle. 'x' and 'y' are
% vectors that contain the x and y-coordinates of the vesicle. 'sigma is a 
% vector that contains the tension on the membrane. 'dsdt' is a vector that 
% contains the derivative of arc length with respect to the 
% parameterization variable. 'dsdt' is computed using Dsdt.m.
%

% Computes (\sigma dx/ds)
DT1=Dtp([x,y]);
x1=sigma.*DT1(:,1)./dsdt;
y1=sigma.*DT1(:,2)./dsdt;

% Differentiates (\sigma dx/ds) with respect to arc length.
DT2=Dtp([x1,y1]);
x2=DT2(:,1)./dsdt;
y2=DT2(:,2)./dsdt;

% Outputs results.
X=real([x2,y2]);
end