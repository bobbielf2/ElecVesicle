function Z=SLP_VV_operator2(f,x,y,SS,CM,FM,mu)


M=length(x);
L=zeros(M,1);
for k=1:M
    L(k)=length(x{k})-1;
end
N=sum(L);
F=zeros(2*N,1);
ls=0;


for k=1:M
    SF{k}=SS{k}*[f{k}(2:end,1);f{k}(2:end,2)];
    SF{k}=[SF{k}([end/2,1:end/2]),SF{k}([end,end/2+1:end])];
end

for k=1:M
    F(ls+1:ls+L(k))=f{k}(2:end,1);
    F(N+ls+1:N+ls+L(k))=f{k}(2:end,2);
    ls=ls+L(k);
end

Z=CM*F;

if ~isempty(FM)
    Z=Z+FM*F;
else
    s.x=zeros(N,1);
    s.spn=zeros(N,1);
    ls=0;
    for k=1:M
        sk.x=x{k}(2:end)+1i*y{k}(2:end);
        sk=quadr(sk);
        s.x(ls+1:ls+L(k))=sk.x;
        s.spn(ls+1:ls+L(k))=sk.sp/L(k);
        ls=ls+L(k);
    end
    t.x=s.x;
    Z=Z+SLPFMM(t,s,F,4);
    slr.x=[s.x+2*pi;s.x-2*pi];
    slr.spn=[s.spn;s.spn];
    Flr=[F(1:end/2);F(1:end/2);F(end/2+1:end);F(end/2+1:end)];
    Z=Z+SLPFMM(t,slr,Flr,4);
end

ls=0;
for k=1:M
    Z(ls+1:ls+L(k))=Z([ls+L(k),ls+1:ls+L(k)-1])+SF{k}(1:end-1,1);
    Z(N+ls+1:N+ls+L(k))=Z([N+ls+L(k),N+ls+1:N+ls+L(k)-1])+SF{k}(1:end-1,2);
    ls=ls+L(k);
end

Z=Z/mu;
end
