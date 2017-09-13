function Z=P_operator(X,xold,yold)
%
% This function is the D operator. 'X' is the input to the operator  with
% the form X=[x1;y1;x2;y2;...], where 'xj' and 'yj' are x and y-coordinates
% of the jth vesicle (not including the last point). 'xold' and 'yold' are
% cells that contain the x and y-coordinates of the vesicles from the
% previous time-step. 'Z' is the output of the D operator with the same
% format as 'X'.
%
% See also:
%   Dsdt, Ds, Dt, VWM
%

% Variable index:
%   M = Number of vesicles.
%   L = A vector that contains the number of discretization points per
%       vesicle (not including the last point).
%   Z = Output of the D operator.
%   Sf = A cell that contains the x and y-coordinates of the vesicles.
%   dsdt = A cell that contains the derivative of arc length with respect
%       to the parameterization variable usually referred to as alpha in
%       literature.
%   dxds = A vector that contains the derivative of x with respect to arc
%       length.
%   dyds = A vector that contains the derivative of y with respect to arc
%       length.
%   Sfs = A vector that contains the derivative of the input x and
%       y-coordinates with respect to arc length.
%   z = A vector that contains the dot product of 'Sfs' with the derivative
%       of the position vector with respect to arc length.
%

% Initializes variables
M=length(xold);
L=zeros(M,1);
Z=zeros(length(X)/2,1);
for k=1:M
    L(k)=length(xold{k})-1;
end
Sf=cell(1,M);
rs=0;
for k=1:M
    Sf{k}=[X(rs+1:L(k)+rs),X(end/2+rs+1:end/2+rs+L(k));X(rs+1),X(end/2+rs+1)];
    rs=L(k)+rs;
end

% Computes the derivative of arc length with respect to the
% parameterization variable.
dsdt=cell(1,M);
for k=1:M
    dsdt{k}=Dsdt(xold{k},yold{k});
end

% Applies the D operator on each vesicle.
rs=0;
for k=1:M
    [dxds,dyds]=Ds(xold{k},yold{k},dsdt{k});
    Sfs=Dtp(Sf{k})./(dsdt{k}*ones(1,2));
    z=dxds.*Sfs(:,1)+dyds.*Sfs(:,2);
    Z(rs+1:rs+L(k))=z(1:L(k));
    rs=rs+L(k);
end