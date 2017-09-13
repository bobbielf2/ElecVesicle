function tau = f_el(s, S, Sprime, E_infty, Vm, q)

beta = 1; %field strength
epsilon = beta*ones(size(s));
[N, M] = size(Vm);

phi_in = zeros(N,M);
phin_in = phi_in;
E_in = phi_in;
E_ex = phi_in;

x = [];
for i=1:length(s)
    x = [x;s{i}.x];
end

%2, compute Dprime
DsVm = zeros(size(Vm(:)));
for i = 1:M
    ind_i = (i-1)*N+(1:N);
    DsVm(ind_i) = diff_fft(Vm(:,i))./s{i}.sp;
end

SDsVm = S*DsVm;

Dprime = zeros(size(SDsVm));
for i = 1:M
    ind_i = (i-1)*N+(1:N);
    Dprime(ind_i) = diff_fft(SDsVm(ind_i))./s{i}.sp;
end

% compute phi & normal deriv of phi on bdry
for i = 1:M
    phi_in(:,i) = -real(conj(E_infty)*s{i}.x) + Vm(:,i)/2 - DLPmatrix(s{i},s{i})*Vm(:,i);
    phin_in(:,i) = -real(conj(E_infty)*s{i}.nx) + q(:,i)/2;
    for j = [1:i-1,i+1:M]
        % phi
        phi_in(:,i) = phi_in(:,i) - lapDevalclose(s{i}.x,s{j},Vm(:,j),'e');
    end
end

phi_in = phi_in + reshape(S*q(:),N,M);
phin_in = phin_in + reshape(Sprime*q(:),N,M) - reshape(Dprime,N,M);

phi_ex = phi_in - Vm;
phin_ex = phin_in - q;

% grad phi = tang deriv * tang + normal deriv * normal
for i = 1:M
    E_in(:,i) = -diff_fft(phi_in(:,i))./s{i}.sp.*s{i}.tang - phin_in(:,i).*s{i}.nx;
    E_ex(:,i) = -diff_fft(phi_ex(:,i))./s{i}.sp.*s{i}.tang - phin_ex(:,i).*s{i}.nx;
end


% compute elec stress

n = zeros(size(E_in)); % bdry normal
for i = 1:M
    n(:,i) = s{i}.nx;
end

% tau = jump of [(bdry normal) dot (maxwell tensor)]
E = E_in;
tau_in = E.*real(conj(E).*n) - abs(E).^2.*n/2;
E = E_ex;
tau_ex = E.*real(conj(E).*n) - abs(E).^2.*n/2;
tau = (tau_in - tau_ex)*diag(epsilon);
tau = [real(tau);imag(tau)];


function s = quadr(s)  % set up periodic trapezoid quadrature on a segment
% Note sign change in normal vs periodicdirpipe.m.  Barnett 4/21/13
N = length(s.x); s.xp = diff_fft(s.x,1); s.xpp = diff_fft(s.x,2);
s.t = (1:N)'/N*2*pi; s.sp = abs(s.xp); s.tang = s.xp./s.sp; s.nx = -1i*s.tang; 
s.cur = -real(conj(s.xpp).*s.nx)./s.sp.^2; s.w = 2*pi/N*s.sp; % speed weights
s.cw = 1i*s.nx.*s.w;  % complex weights (incl complex speed)
s.a = mean(s.x);

function [A An] = DLPmatrix(t,s,a) % double-layer kernel matrix & targ n-deriv
% t = target seg (x,nx cols), s = src seg, a = optional translation of src seg
% No jump included on self-interaction (ie principal-value integral)
if nargin<3, a = 0; end
N = numel(s.x); M = numel(t.x);
d = repmat(t.x, [1 N]) - repmat(s.x.' + a, [M 1]);    % C-# displacements mat
ny = repmat(s.nx.', [M 1]);      % identical rows given by src normals
A = (1/2/pi) * real(ny./d);      % complex form of dipole
if numel(s.x)==numel(t.x) & max(abs(s.x+a-t.x))<1e-14
  A(diagind(A)) = -s.cur/4/pi; end           % self? diagonal term for Laplace
A = A .* repmat(s.w(:)', [M 1]);
if nargout>1 % deriv of double-layer. Not correct for self-interaction.
  csry = conj(ny).*d;              % (cos phi + i sin phi).r
  nx = repmat(t.nx, [1 N]);        % identical cols given by target normals
  csrx = conj(nx).*d;              % (cos th + i sin th).r
  r = abs(d);    % dist matrix R^{MxN}
  An = -real(csry.*csrx)./(r.^2.^2)/(2*pi);
  An = An .* repmat(s.w(:)', [M 1]);
end

function [A An] = SLPmatrix(t,s,a) % single-layer kernel matrix & targ n-deriv
% t = target seg (x,nx cols), s = src seg, a = optional translation of src seg
% No jump included on self-interaction of derivative (ie PV integral).
if nargin<3, a = 0; end
N = numel(s.x); M = numel(t.x);
d = repmat(t.x, [1 N]) - repmat(s.x.' + a, [M 1]);    % C-# displacements mat
A = -(1/2/pi) * log(abs(d)) .* repmat(s.w(:)', [M 1]);  % infty for self-int
if nargout==2                      % apply D^T
  nx = repmat(-t.nx, [1 N]);       % identical cols given by -targ normals
  An = (1/2/pi) * real(nx./d);     % complex form of dipole. Really A1 is An
  if numel(s.x)==numel(t.x) & max(abs(s.x+a-t.x))<1e-14
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
