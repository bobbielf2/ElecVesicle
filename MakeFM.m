function FM=MakeFM(s)

M=length(s);
L=zeros(M,1);
for k=1:M
   L(k)=size(s{k},1);
end
FM=zeros(2*sum(L));
ss=cell(1,M);
for k=1:M
ss{k}.x=s{k}(:,1)+1i*s{k}(:,2);
ss{k}=quadr(ss{k});
end

lt=0;
for j=1:M
    ls=0;
    for k=1:M
        SLP=SLPmatrix(ss{j},ss{k},1);
        if j==k
            SLP(logical([eye(size(SLP)/2),eye(size(SLP)/2);eye(size(SLP)/2),eye(size(SLP)/2)]))=0;
        end
        FM(lt+1:lt+L(j),ls+1:ls+L(k))=SLP(1:end/2,1:end/2);
        FM(end/2+lt+1:end/2+lt+L(j),ls+1:ls+L(k))=SLP(end/2+1:end,1:end/2);
        FM(lt+1:lt+L(j),end/2+ls+1:end/2+ls+L(k))=SLP(1:end/2,end/2+1:end);
        FM(end/2+lt+1:end/2+lt+L(j),end/2+ls+1:end/2+ls+L(k))=SLP(end/2+1:end,end/2+1:end);
        ls=ls+L(k);
    end
    lt=lt+L(j);
end