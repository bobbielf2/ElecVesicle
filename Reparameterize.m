function [X,Y,SIG]=Reparameterize(X,Y,SIG,k)
%
% This function reparameterizes the vesicles so that discretization points
% are approximately equally spaced. 'X', 'Y', 'SIG', and 'S' are the same
% as in VWM.m. 'k' is the current time-step within VWM.m.
%
% See also:
%   Dt, INTt, It, VWM
%

% Variable index:
%   M = Number of vesicles.
%   n = The number of discretization points for the jth vesicle.
%   t = Vesicle parameterization variable. This is usually alpha in
%       literature.
%   ds = The derivative of arc length with respect to alpha for the jth
%       vesicle.
%   s = Arc length of the jth vesicle and also the desired arc length
%       values for the new discretization.
%   st = Integral of ds with respect to alpha.
%   ti = Alpha values for new discretization points.
%

M=length(Y{1});
for j=1:M
    n=length(X{k}{j});
    t=linspace(0,2*pi,n)';
    ds=sqrt(Dtp(X{k}{j}).^2+Dtp(Y{k}{j}).^2);
    s=sum(ds(1:n-1))*2*pi/(n-1);
    s=linspace(0,s,n)';
    st=INTt(ds);
    
%     tic
    tj=zeros(n,1);
    tj(1)=0;
    tj(n)=2*pi;
    jj=1;
    for ii=2:n-1
       while ii<n
          if s(ii)>=st(jj)&&s(ii)<st(jj+1)
              tj(ii)=2*pi/(n-1)*(jj-1+(s(ii)-st(jj))/(st(jj+1)-st(jj)));
              break
          else
              jj=jj+1;
          end
       end
    end
    ti=tj;
%     toc
%     
%     tic
%     % Computes new alpha values
%     ST=repmat(s',n,1)-repmat(st,1,n);
%     STp=ST;
%     STp(STp<0)=inf;
%     STn=ST;
%     STn(STn>=0)=-inf;
%     [minv, minr]=min(STp,[],1);
%     [maxv, maxr]=max(STn,[],1);
%     ti=(s-st(minr))./(st(maxr)-st(minr)).*(t(maxr)-t(minr))+t(minr);
%     ti(1)=0;
%     ti(n)=2*pi;
%     toc
    
    % Interpolates 'X', 'Y', and 'SIG' at these values
    X{k}{j}=Itp(X{k}{j},ti);
    Y{k}{j}=Itp(Y{k}{j},ti);
    if k>1
        SIG{k-1}{j}=Itp(SIG{k-1}{j},ti);
    end
end
end