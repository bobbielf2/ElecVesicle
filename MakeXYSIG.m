function [X,Y,SIG]=MakeXYSIG(filename)

load(filename)
X=cell(1,k+1);
Y=cell(1,k+1);
SIG=cell(1,k);
X{1}=X_1;
Y{1}=Y_1;
for k=1:k
    eval(['X{k+1}=X_',num2str(k+1),';'])
    eval(['Y{k+1}=Y_',num2str(k+1),';'])
    eval(['SIG{k}=SIG_',num2str(k),';'])
end