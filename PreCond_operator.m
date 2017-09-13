function Z=PreCond_operator(alpha,beta,gamma,PreInv,M,LS)

n=sum(LS);
Z=zeros(3*n,1);
ls=0;
for k=1:M
    Z([ls+1:ls+LS(k),end/3+ls+1:end/3+ls+LS(k),2*end/3+ls+1:2*end/3+ls+LS(k)])=PreInv{k}*[alpha(ls+1:ls+LS(k));beta(ls+1:ls+LS(k));gamma(ls+1:ls+LS(k))];
    ls=ls+LS(k);
end