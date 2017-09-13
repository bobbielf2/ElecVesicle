function F=Fsig_operator(gamma,x,y)

M=length(x);
L=zeros(M,1);
for k=1:M
    L(k)=length(x{k})-1;
end
if isequal(class(gamma),'double')
   Gamma=cell(1,M);
   ls=0;
   for k=1:M
       Gamma{k}=[gamma(ls+1:ls+L(k));gamma(ls+1)];
       ls=ls+L(k);
   end
else
   Gamma=gamma;
end

F=cell(1,M);
for k=1:M
   dsdt=Dsdt(x{k},y{k});
   F{k}=fsig(x{k},y{k},Gamma{k},dsdt);
end