function L = L_operator_Vm(Sprime,sig_i,sig_e)
% Input: Sprime = normal deriv of laplace SLP

% L operator
M = length(sig_i);
N = length(Sprime)/M;
eta = (sig_i-sig_e)./(sig_i+sig_e);
eta = eta(:);
Eta = repmat(kron(eta,ones(N,1)),1,M*N); % Eta matrix

L = eye(size(Sprime))/2+Eta.*Sprime;



