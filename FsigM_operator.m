function Z=FsigM_operator(x,y)

L=length(x)-1;
dsdt=Dsdt(x,y);
EYE=eye(L);
EYE=EYE([1:end,1],:);

[dxds,dyds]=Ds(x,y,dsdt);
DS=Dtp(EYE)./(dsdt*ones(1,L));
DS=DS(1:end-1,:);
Z=[(ones(L,1)*dxds(1:end-1).').*DS;(ones(L,1)*dyds(1:end-1).').*DS];