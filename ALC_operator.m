function IE=ALC_operator(x,y,ArcL,dt)

M=length(x);
L=zeros(M,1);
for k=1:M
    L(k)=length(x{k})-1;
end
IE=zeros(sum(L),1);

% This section provides a correction to the arc length in order to
% insure that inextensibility is maintained for long-time
% simulations. (This section of code is optional and may be
% commented out.)
ls=0;
for j=1:M
    IE(ls+1:ls+L(j))=(ArcL(j)-ArcLength(x{j},y{j},2*pi))/ArcLength(x{j},y{j},2*pi)/dt;
    ls=ls+L(j);
end