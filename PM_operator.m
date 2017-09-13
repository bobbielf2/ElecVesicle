function Z=PM_operator(x,y)
%
% This function is the D operator. 'X' is the input to the operator  with
% the form X=[x1;y1;x2;y2;...], where 'xj' and 'yj' are x and y-coordinates
% of the jth vesicle (not including the last point). 'xold' and 'yold' are
% cells that contain the x and y-coordinates of the vesicles from the
% previous time-step. 'Z' is the output of the D operator with the same
% format as 'X'.
%

% Initializes variables
L=length(x)-1;

% Computes the derivative of arc length with respect to the
% parameterization variable.
dsdt=Dsdt(x,y);

% Applies the D operator on each vesicle.
EYE=eye(L);
EYE=EYE([1:end,1],:);
[dxds,dyds]=Ds(x,y,dsdt);
Sfs=Dtp(EYE)./(dsdt*ones(1,L));
Z=[(dxds*ones(1,L)).*Sfs,(dyds*ones(1,L)).*Sfs];
Z=Z(1:end-1,:);