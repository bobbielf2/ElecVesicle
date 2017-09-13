function [vc vcp] = cauchycompeval(x,s,vb,side)
% CAUCHYCOMPEVAL - global compensated evaluation int/exterior Cauchy integral
%
% v = cauchycompevalext(x,s,vb,side) approximates the holomorphic
%  function v at targets x exterior to a closed curve s on which its boundary
%  values vb (v^+ for exterior case, v^- for interior) are given at nodes.
%  "side" specifies if the targets in interior or exterior to the curve.
%  For the exterior case, uniqueness of the holomorphic function is insured
%  by assuming the limit as |z|->infty is zero.
%
% [v vp] = cauchycompevalext(x,s,vb,side) also returns v's complex derivative
%
% For both values and derivatives, true barycentric forms are used that ensure
% machine accuracy for arbitrarily small (or zero) target-node distances.
%
% Inputs:
% x = row or col vec of M targets in complex plane
% s = closed curve quadrature struct containing N nodes s.x, and s.w
%     "speed weights", s.nx unit normals, s.a interior pt far from the bdry.
% vb = row or col vec of N boundary values of holomorphic function v
% side = 'i' or 'e' appropriate for if targets interior or exterior to curve.
% Outputs:
% v  = row vec of homolorphic function v at the M targets
% vp = (optional) row vec of complex first derivative v' at the M targets
%
% Notes; 1) algorithm is that of Ioakimidis et al BIT 1991 for interior, and,
% for the exterior, a modified version using 1/(z-a) in place of the function 1.
% For the derivative, the formula is mathematically the derivative of the
% baycentric formula of Schneider-Werner 1986 as described in Berrut et al 2005,
% however, since this derivative is not a true barycentric formula, a hack is
% needed to compute the difference v_j-v(x) in a form where roundoff error
% cancels correctly. The exterior case is slightly more intricate. When
% a target coincides with a node j, the true value v_j (or Schneider-Werner
% formula for v'(x_j)) is used.
% 2) In order to vectorize in both the node and target indices, we use a lot
% of RAM, O(NM), which is not really necessary. For now, the user should group
% their targets into smaller blocks done with separate callls if the RAM use
% is too high.
%
% Alex Barnett 10/22/13 based on cauchycompevalint, pole code in lapDevalclose.m
% 10/23/13 node-targ coincidences fixed, exterior true barycentric discovered.

mindist = 1e-3;    % dist param where $O(N^2.M) deriv bary form switched on.
                   % Roughly deriv errors are then limited to emach/mindist
                   % The only reason to decrease mindist is if lots of v close
                   % nodes on the curve with lots of close target points.
cw = 1i*s.nx.*s.w;                    % complex speed weights
N = numel(s.x); M = numel(x);

if nargout==1  % no deriv wanted... (note sum along 1-axis faster than 2-axis)
  comp = repmat(cw(:), [1 M]) ./ (repmat(s.x(:),[1 M]) - repmat(x(:).',[N 1]));
  if side=='e', pcomp = comp .* repmat(1./(s.x(:)-s.a), [1 M]);
  else pcomp = comp; end  % pcomp are weights and bary poles appearing in J0
  I0 = sum(repmat(vb(:),[1 M]).*comp); J0 = sum(pcomp); % Ioakimidis notation
  vc = I0./J0;                        % bary form
  if side=='e', vc = vc./(x(:).'-s.a); end         % correct w/ pole
  [jj ii] = ind2sub(size(comp),find(~isfinite(comp))); % node-targ coincidences
  for l=1:numel(jj), vc(ii(l)) = vb(jj(l)); end % replace each hit w/ corresp vb
  
else           % 1st deriv wanted...
  invd = 1./(repmat(s.x(:),[1 M]) - repmat(x(:).',[N 1])); % 1/displacement mat
  comp = repmat(cw(:), [1 M]) .* invd;
  if side=='e', pcomp = comp .* repmat(1./(s.x(:)-s.a), [1 M]);
  else pcomp = comp; end  % pcomp are weights and bary poles appearing in J0
  I0 = sum(repmat(vb(:),[1 M]).*comp); J0 = sum(pcomp);
  if side=='e', prefac = 1./(x(:).'- s.a); else prefac = 1; end
  vc = prefac .* I0./J0;                        % bary form (poss overall pole)
  dv = repmat(vb(:),[1 M]) - repmat(vc(:).',[N 1]); % v value diff mat
  [jj ii] = ind2sub(size(invd),find(abs(invd) > 1/mindist)); % bad pairs indices
  if side=='e' % exterior:
    for l=1:numel(jj), j=jj(l); i=ii(l); % loop over node-targ pairs too close
      p = sum(comp(:,i).*(vb(j)./(s.x(:)-s.a)-vb(:)./(s.x(j)-s.a))) / sum(pcomp(:,i)); % p is v_j - (x-a)/(yj-a)*v(x) for x=i'th target
      dv(j,i) = prefac(i) * ((s.x(j)-s.a)*p + (x(i)-s.x(j))*vb(j));
    end % pole-corrected for dv, gives bary stability for close node-targ pairs
  else  % interior:
    for l=1:numel(jj), j=jj(l); i=ii(l); % loop over node-targ pairs too close
      dv(j,i) = sum(comp(:,i).*(vb(j)-vb(:))) / sum(pcomp(:,i));
    end % bary for dv, gives bary stability for close node-targ pairs
  end
  vcp = prefac .* sum(dv.*comp.*invd) ./ J0; % bary form for deriv
  [jj ii] = ind2sub(size(comp),find(~isfinite(comp))); % node-targ coincidences
  for l=1:numel(jj), j=jj(l); i=ii(l); % loop over hitting node-targ pairs
    vc(i) = vb(j);                     % replace each hit w/ corresp vb
    notj = [1:j-1, j+1:N];             % Schneider-Werner form for deriv @ node:
    vcp(i) = -sum(cw(notj).*(vb(j)-vb(notj))./(s.x(j)-s.x(notj)))/cw(j);
  end
end
