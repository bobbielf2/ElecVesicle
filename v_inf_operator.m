function v_inf = v_inf_operator(alpha,beta,x,shear)

M=length(x);
L=zeros(M,1);
for k=1:M
    L(k)=length(x{k})-1;
end
N = sum(L);

v_inf = zeros(2*N,1);

if shear == 0, return, end

ls = 0;
for k=1:M
    v_inf(ls+1:ls+L(k)) = shear*beta([ls+2:ls+L(k),ls+1]);
    v_inf(N+ls+1:N+ls+L(k))=0;
    ls=ls+L(k);
end




