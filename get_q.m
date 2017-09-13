function q = get_q(s,Vm,S,L,E_infty,sigma_i,sigma_e,eps)
% compute the jump of the normal derivative of the electric potential
% Input:
%  s = cell array, each cell is a structure containing quadrature
%      information of a single vesicle. e.g. s{1}.x = cos(t) + 1i*sin(t);
%      s{1} = quadr(s{1});

% eps = 1e-9;
[N,M] = size(Vm);
Vm = Vm(:);

%buid rhs vector
eta = (sigma_i - sigma_e)./(sigma_i + sigma_e); 
Eta = kron(diag(eta),eye(N)); % Eta matrix
Nx = []; % outward normals on bdry
for i = 1:M, Nx = [Nx;s{i}.nx]; end

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

rhs = Eta*(real(conj(E_infty)*Nx) + Dprime);

q = gmres(L,rhs,[],eps,length(rhs));
q = reshape(q,[],M);



function [A, An] = SLPmatrix(t,s,a) % single-layer kernel matrix & targ n-deriv
% t = target seg (x,nx cols), s = src seg, a = optional translation of src seg
% No jump included on self-interaction of derivative (ie PV integral).
if nargin<3, a = 0; end
N = numel(s.x); M = numel(t.x);
d = repmat(t.x, [1 N]) - repmat(s.x.' + a, [M 1]);    % C-# displacements mat
A = -(1/2/pi) * log(abs(d)) .* repmat(s.w(:)', [M 1]);  % infty for self-int
if nargout==2                      % apply D^T
  nx = repmat(-t.nx, [1 N]);       % identical cols given by -targ normals
  An = (1/2/pi) * real(nx./d);     % complex form of dipole. Really A1 is An
  if numel(s.x)==numel(t.x) && max(abs(s.x+a-t.x))<1e-14
    An(diagind(An)) = -s.cur/4/pi; end  % self? diagonal term for Laplace
  An = An .* repmat(s.w(:)', [M 1]);
end

    
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
