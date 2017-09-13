function PreInv=PreInv_operator(x,y,SS,kb,dt,mu)

M=length(x);
PreInv=cell(1,M);
for k=1:M
    L=length(x{k})-1;
    PM=PM_operator(x{k},y{k});
    FbM=FbM_operator(x{k},y{k},kb);
    FbM=FbM([2:end/2,1,end/2+2:end,end/2+1],:);
    FsigM=FsigM_operator(x{k},y{k});
    FsigM=FsigM([2:end/2,1,end/2+2:end,end/2+1],:);
    SS{k}=SS{k}([end/2,1:end/2-1,end,end/2+1:end-1],:)/mu;
    PreInv{k}=inv([PM,zeros(L);[eye(2*L),zeros(2*L,L)]-SS{k}*[dt*FbM,FsigM]]);
end