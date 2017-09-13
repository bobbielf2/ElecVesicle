function Z = LinOP_Vm(Vm,s,S,L,Cm,Gm,dt,sig_i,sig_e)
% Input: S = laplace SLP mtx
%        Sprime = normal deriv of S
N = length(s{1}.x);
M = length(s);
lambda = sig_i*sig_e./(sig_i+sig_e);

%1, compute L's
CM = kron(Cm(:),ones(N,1)); % Cm matrix
GM = kron(Gm(:),ones(N,1)); % Gm matrix
Lambda = kron(lambda(:),ones(N,1)); % Lambda matrix

Z = (CM+dt*GM).*(L*Vm);

%2, compute Dprime
DsVm = zeros(size(Vm));
for i = 1:M
    ind_i = (i-1)*N+(1:N);
    DsVm(ind_i) = diff_fft(Vm(ind_i))./s{i}.sp;
end

SDsVm = S*DsVm;

Dprime = zeros(size(Vm));
for i = 1:M
    ind_i = (i-1)*N+(1:N);
    Dprime(ind_i) = diff_fft(SDsVm(ind_i))./s{i}.sp;
end

Z = Z - dt*Lambda.*Dprime;


    
function df = diff_fft(f,m)
if nargin<2, m = 1; end
N = size(f,1);
if(mod(N,2)==1),warning('DmFT:size',['Are you sure about using FFT ' ...
        'for differentiation of a vector \n with odd ' ...
        'length?'])
end
k = 1i*[(0:N/2-1)  0 (-N/2+1:-1)]';

df = f;

for j = 1:m
    df = ifft(k.*fft(df));
end

