function y=Fltr(x)
%
% This function filters out the 1/3 highest frequencies of 'x'. The first
% point of 'x' must be in the in the first and last row of 'x'.
%

[m,p]=size(x);
alpha=linspace(0,1,m)'*(x(m,:)-x(1,:));
m=m-1;
X=fft(x(1:m,:)-alpha(1:m,:));
Y=zeros(m,p);
n=floor((floor(2*m/3)-1)/2);
Y(1:n+1,:)=X(1:n+1,:);
Y(m-n+1:m,:)=X(m-n+1:m,:);
y=ifft(Y)+alpha(1:m);
y=real([y;y(1,:)+alpha(m+1,:)]);