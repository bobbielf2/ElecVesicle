function S=SLP_self_operator(x,y)

M=length(x);
S=cell(1,M);
for k=1:M
   S{k}=SLPselfmatrix([x{k},y{k}]); 
end