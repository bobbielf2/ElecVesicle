function [X,Y]=Filter(X,Y,k)

M=length(Y{1});
for j=1:M
    X{k}{j}=Fltr(X{k}{j});
    Y{k}{j}=Fltr(Y{k}{j});
end