function [Z, Zp] = SLP_VV_lap(s)
% compute the laplace SLP for vesicle-vesicle interactions

%Step 1, find h, will perform closeeval for vesicle distance<5h
M = length(s);
N = length(s{1}.x);
[D,h] = Findh(s); %min dist b/w nodes on a vesicle

%Step 2, construct SLP matrix
Z = zeros(M*N);
Zp = Z;

for i = 1:M
    ind_i = (i-1)*N+1:i*N;
    for j = 1:M
        ind_j = (j-1)*N+1:j*N;
        %construct block S_ij, the interaction between vesicle i and j
        if i==j
            %self matrix
            Z(ind_i,ind_j) = kernelSLapk(s{i});
            [~, Zp(ind_i,ind_j)] = SLPmatrix(s{i},s{i});
            
        elseif D(i,j)<5*h %distant between vesicle i and j < 5h
            %close evaluation
            [Z(ind_i,ind_j),Zpx,Zpy] = lapSevalcloseF(s{i}.x, s{j}, [], 'e');
            Zp(ind_i,ind_j) = repmat(real(s{i}.nx),1,N).*Zpx+repmat(imag(s{i}.nx),1,N).*Zpy;
        else
            %Nystrom matrix
            [Z(ind_i,ind_j),Zpx,Zpy] = lapSevalcloseF(s{i}.x, s{j}, [], 'e');
            Zp(ind_i,ind_j) = repmat(real(s{i}.nx),1,N).*Zpx+repmat(imag(s{i}.nx),1,N).*Zpy;
        end
    end
end



function [D,h] = Findh(s)
M = length(s);
N = length(s{1}.x);
D = zeros(M);
h=Inf;
for i=1:M
    X = s{i}.x([1:end,1]);
    d = ArcLength(real(X),imag(X),2*pi)/N;
    h = min(h,d);
    for j = 1:M
        if i==j
            D(i,j) = 0;
        else
            dist = abs(repmat(s{i}.x,1,N) - repmat(s{j}.x.',N,1));
            D(i,j) = min(dist(:));
        end
    end
end

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
  if numel(s.x)==numel(t.x) & max(abs(s.x+a-t.x))<1e-14
    An(diagind(An)) = -s.cur/4/pi; end  % self? diagonal term for Laplace
  An = An .* repmat(s.w(:)', [M 1]);
end