function f=Fdt_operator(alpha,beta,gamma,Fb,Fsig,dt)

fsig=Fsig(gamma);
fb=Fb(alpha,beta);
M=length(fsig);
f=cell(1,M);
for k=1:M
   f{k}=fsig{k}+dt*fb{k}; 
end