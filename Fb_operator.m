function F=Fb_operator(alpha,beta,x,y,kb)

M=length(x);
L=zeros(M,1);
for k=1:M
    L(k)=length(x{k})-1;
end
if isequal(class(alpha),'double')
   Alpha=cell(1,M);
   Beta=Alpha;
   ls=0;
   for k=1:M
       Alpha{k}=[alpha(ls+1:ls+L(k));alpha(ls+1)];
       Beta{k}=[beta(ls+1:ls+L(k));beta(ls+1)];
       ls=ls+L(k);
   end
else
   Alpha=alpha;
   Beta=beta;
end

F=cell(1,M);
for k=1:M
   dsdt=Dsdt(x{k},y{k});
   F{k}=fb(Alpha{k},Beta{k},dsdt,kb);
end