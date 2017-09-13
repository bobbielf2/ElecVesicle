function v_inf = V_inf(X,Y,shear)
M = length(X);
for k=1:M
    L(k)=length(X{k})-1;
end
N = sum(L);
v_inf=zeros(sum(L)*2,1);

if shear == 0, return, end

rs = 0;
for k=1:M
        v_inf(rs+1:rs+L(k))=shear*Y{k}(2:end);
        v_inf(rs+N+1:rs+N+L(k))=0;
        rs=rs+L(k);

end

end
