function S =  kernelSLapk(s)
% Laplace SLP self operator via Kress quad


N = length(s.x);
d = repmat(s.x, [1 N]) - repmat(s.x.', [N 1]);
r = sqrt(d.*conj(d));    % dist matrix R^{MxN}
logpart =  -(1/2/pi) * log(r);
S1 = -1/4/pi;
logpart = logpart  - S1.*circulant(log(4*sin(pi*(0:N-1)/N).^2));
logpart(diagind(logpart)) = -log(s.sp)/2/pi;   % diag vals propto curvature
m = 1:N/2-1; Rjn = -2*pi*ifft([0 1./m 2/N 1./m(end:-1:1)]);
logpart = (circulant(Rjn).*S1 + logpart*(2*pi/N)) .* repmat(s.sp.', [N 1]);
% logpart = logpart / 2;  % hack the prefactor

% S = logpart/4/pi; % WHY is the factor 4pi?
S = logpart;