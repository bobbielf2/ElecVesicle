function avgS = avgStressElec(X,u,sig,f_el,prams)
  
  if(nargin==0), avgS = avgStressTest(); return; end
%   avgS = zeros(2); area=0;  
  %area = calcArea(X(:));
  
  n = size(X,1)/2;
  [kap, tang, sa] = curveProp(X);
  tang = reshape(tang,[],2);
  nor = [tang(:,2) -tang(:,1)];
  SA = reshape(2*pi/n*repmat(sa',4,1),2,2,n);
  X = reshape(X,[],2);
  
%   cn = -prams.kappa*kap.^2; ct = -sig;
%   
%   val1 = tensorProd([cn cn].*nor,nor) + tensorProd([ct ct].*tang,tang)- tensorProd(f_el,X);

  sa = sqrt(DmFT(X(:,1),[]).^2 + DmFT(X(:,2),[]).^2);
  x1 = DmFT(X,1,[],1./sa);
  f = -prams.kappa*DmFT(x1,3)+ DmFT([sig.*x1(:,1), sig.*x1(:,2)], 1) - f_el;
  
  val3 = tensorProd(f,X);
%   disp(sum((val1-val3).*SA,3));
  
  %uInf = prams.vInf(X(:)); u = u-uInf;
  u = reshape(u,[],2);
  
  val2 = tensorProd(u,nor)+tensorProd(nor,u); 
  
  val = (-val3+(prams.viscCont-1)*val2).*SA;  %-val1 because of the signs.
  avgS = sum(val,3);
  % This is the relevant term for the effective viscosity
  avgS = avgS(1,2);
  
function prod = tensorProd(x,y)
% The input format is one/two column vectors 
  if(size(x,2)==1)
    x = reshape(x,[],2);
    y = reshape(y,[],2);
  elseif(size(x,2)~=2)
    disp('check sizes');
    prod = [];
    return;
  end
  
  prod(1,1,:) = x(:,1).*y(:,1);
  prod(1,2,:) = x(:,1).*y(:,2);
  prod(2,1,:) = x(:,2).*y(:,1);
  prod(2,2,:) = x(:,2).*y(:,2);
  
function df = DmFT(f, m, kIn, saIn)
%DmFT(f,m) Computes the mth derivative of a periodic function f using
%fourier transform. By calling DmFT(f,m,K,Sa) one can overrride the
%differentiation coefficients with the desired ones by passing extra
%arguments K and the jacobian Sa.
%
%EXAMPLE:
%  n = 64;
%  X = boundary(n,'nv',1); X = reshape(X,n,2);
%  DX  = DmFT(X,1);
%  DDX = DmFT(X,2);
%
%  sa = sqrt(dot(DX,DX,2));
%  tang = DX./[sa sa];
%  kap = (DX(:,1).*DDX(:,2) - DX(:,2).*DDX(:,1))./sa.^3;
% 
%  viewer(X(:),tang(:)); 
%

  persistent k sa

  if(nargin == 1 || isempty(m)), m=1 ;end
  if(nargin >= 3), k = kIn; end
  if(nargin == 4), sa = saIn; end

  if(isempty(k))
    N = size(f,1);
    if(mod(N,2)==1),warning('DmFT:size',['Are you sure about using FFT ' ...
                          'for differentiation of a vector \n with odd ' ...
                          'length?'])
    end
    k = 1i*[(0:N/2-1)  0 (-N/2+1:-1)]';
  end

  if(isempty(sa)), sa = ones(size(f,1),1);end

  df = f;

  if(~isempty(f))
    [n, col] = size(f);
    
    if(size(sa,2)~=col)
      S = repmat(sa,1,col);
    else
      S = sa;  
    end
    
    K = repmat(k,1,col);
    
    for j = 1:m
      df = S.*real(ifft(K.*fft(df)));
    end
  end
  
function [kap, tang, sa] = curveProp(X)
% CURVEPROP(X) returns the curvature, tangent vector and the Jacobian of
% each column of the matrix X. Each column of X should be a closed curve
% defined in plane. 
%
% EXAMPLE:
%    n = 128; nv = 3;
%    X = boundary(n,'nv',nv,'curly');
%    [k t s] = curveProp(X);
%    viewer(X,t);
%    [theta,rho] = cart2pol([X(1:n,2);X(1,2)],[X(n+1:2*n,2);X(n+1,2)]);
%    figure; polar(theta,[s(:,2);s(1,2)],'r'); hold on; polar(theta,rho,'b'); 
%    legend('The Jacobian','Vesicle boundary','Location','Best');
%    title('The Jacobian of boundary points of the vesicle in polar coordinate')
%
[n, nv] = size(X); n = n/2; 

Y = X(1+n:2*n,:); X = X(1:n,:); 
DX(:,1) = reshape(DmFT(X,1,[],[]),[],1);
DX(:,2) = reshape(DmFT(Y,1,[],[]),[],1);

DDX(:,1) = reshape(DmFT(X,2,[],[]),[],1);
DDX(:,2) = reshape(DmFT(Y,2,[],[]),[],1);

sa = sqrt(dot(DX,DX,2)); sa(sa==0) = 1;
t = DX./[sa sa];
tang(1:n,:) = reshape(t(:,1),n,nv);
tang(n+1:2*n,:) = reshape(t(:,2),n,nv);
kap = (DX(:,1).*DDX(:,2) - DX(:,2).*DDX(:,1))./sa.^3;
kap = reshape(kap,n,nv);
sa = reshape(sa,n,nv);
  
  
function avgS = avgStressTest()
  
  n = 64;
  gg = (1:n)'*2*pi/n;
  X = [cos(gg);sin(gg)];
  sig = ones(n,1);
  f_el = [-sin(gg),cos(gg)];
  u = zeros(2*n,1);
  prams.kappa = 0;
  prams.viscCont = 1;
  
  avgS = avgStressElec(X,u,sig,f_el,prams);
%   [avgS area] = avgStress(X,u,sig,prams);
  
%   n = 64;
%   gg = (1:n)'*2*pi/n;
%   xv = [cos(gg);sin(gg)];
%  
%   tJ = xv; [avgS area] = avgStress([],[],[],xv,tJ);
%   disp(['circle with f = x             : error = ' num2str(norm(avgS-pi*eye(2)))]);
%   
%   tJ(n+1:end) = 0; [avgS area] = avgStress([],[],[],xv,tJ);
%   disp(['circle with f_x = x, f_y = 0  : error = ' num2str(norm(avgS-[pi 0;0 0]))]);
%   
%   tJ(1:n) = -xv(n+1:end); tJ(n+1:2*n) = xv(1:n); 
%   [avgS area] = avgStress([],[],[],xv,tJ);
%   disp(['circle with f_x = -y, f_y = x : error = ' num2str(norm(avgS-[0 -pi;pi 0]))]);
%   
%   fRand = rand(2,1);
%   tJ = zeros(size(xv)); tJ([1 n+1]) = fRand;
%   [avgS area] = avgStress([],[],[],xv,tJ);
%   refS = 2*pi/n*fRand*xv([1 n+1])';
%   disp('circle with f_x = 0, f_y = 0,'); 
%   disp(['and a random nonzero vector   : error = ' num2str(norm(avgS- refS))]);
%   
%   
%   xv = [-2+cos(gg);sin(gg)];[3+2*cos(gg);1+2*sin(gg)];
%   fRand = rand(2,1);
%   tJ = zeros(size(xv)); tJ([1 n+1]) = fRand;
%   [avgS area] = avgStress([],[],[],xv,tJ);
%   refS = 2*pi/n*fRand*xv([1 n+1])';
%   disp(['two circles, same as above    : error = ' num2str(norm(avgS-refS))]);  
%   

% function [avgS area] = avgStress(M,domain,mu,Xv,tractionJump)
% 
%   if(nargin==0), avgStressTest(); return; end
%   
%   avgS = zeros(2);
%   area = [];
% 
%   if(~isempty(domain))
%     surfaceTrac = evalTraction(M,domain,mu);
%     pIndex(1) = 0;
%     for i=1:length(M)
%       pIndex(i+1) = pIndex(i)+M(i);
%     end
%     N = size(domain,2)-1;
% 
%     for ind = 1:N+1
%       Sn = surfaceTrac(pIndex(ind)+1:pIndex(ind+1),:);
%       X = domain(ind).X;
%       sa = domain(ind).jacob;
% 
%       T = domain(ind).h*[Sn(:,1).*X(:,1).*sa Sn(:,2).*X(:,1).*sa ...
%                          Sn(:,1).*X(:,2).*sa Sn(:,2).*X(:,2).*sa];
% 
%       avgS = avgS + reshape(sum(T),2,2);
%     end
%     area = calcArea(domain);
%   end
% 
%   if(~isempty(Xv))
%     n = size(Xv,1)/2;
%     [trash1 trash2 saAll] = curveProp(Xv);
%     for ind = 1:size(Xv,2)
%       tj = reshape(tractionJump(:,ind),[],2);
%       X = reshape(Xv(:,ind),[],2);
%       sa = saAll(:,ind);
% 
%       T = 2*pi/n*[tj(:,1).*X(:,1).*sa tj(:,2).*X(:,1).*sa ...
%                   tj(:,1).*X(:,2).*sa tj(:,2).*X(:,2).*sa];
% 
%       avgS = avgS + reshape(sum(T),2,2);
%     end
%     
%     if(isempty(area)), area = calcArea(Xv);end
%   end
