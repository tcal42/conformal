function [f, finv] = conformal(bndry, ctr, varargin)
%CONFORMAL  Chebfun conformal mapping of a smooth domain to the unit disk.
%
%  [F, FINV] = conformal(BNDRY, CTR) returns function handles F and FINV
%  for the forward and inverse maps.  Here BNDRY is a smooth periodic
%  complex chebfun on [0, 2*PI] (a "trigfun") defining the boundary of the
%  domain, and CTR is the point in the domain to be mapped to 0.
%
%  [F, FINV] = conformal(...,'plot') plots the result too.
%
%  Example:
%  
%  circle = chebfun('exp(1i*t)',[0 2*pi],'trig');
%  bndry = real(circle) + .6i*imag(circle); ctr = 0.5;
%  [f, finv] = conformal(bndry, ctr, 'plot');
%
%  The method used is based on the Kerzman-Stein integral equation, with
%  the output maps represented by AAA rational approximants.  This code was
%  written in 2017 by Trevor Caldwell, Anne Greenbaum, and Nick Trefethen.

%% Construct trigfun for the boundary correspondence function.  
% bcfun maps the arc length along the boundary of Omega (scaled to [0, 2*pi]) 
% to the conformal image on the unit circle.

%%
% Invert the bndry function once and for all
dbndry = diff(bndry);
s = cumsum(abs(dbndry));
S = s(end);
dbndry_arg = angle(dbndry);
dbndry_arg = unwrap(dbndry_arg-dbndry_arg(0));
clockwise = (diff(dbndry_arg([0 end])) < 0);
% This is where the time gets spent! :
if clockwise
  u = inv(S-s); 
else
  u = inv(s);
end
minu = u(0); maxu = u(end);
bndry2 = newDomain(bndry,[minu maxu]);
W = bndry2(u);
Wprime = diff(W);

if nargin>2
  clf
  subplot(1,2,1)
  plot(bndry,'k'), hold on
  axis square, axis(1.3*[-1 1 -1 1])
  set(gca,'fontsize',8), drawnow
end

bcfun = chebfun(@(t) bcfun_eval(t, bndry2, W, Wprime, S, ctr), [0 2*pi], 'trig', ...
        'sampleTest', 0, 'resampling', 'on'); 
bcfun = chebfun(bcfun.values,[0 2*pi],'trig');
param_pts = pi*(1+bcfun.points); % equispaced pts in [0,2*pi)
domain_pts = bndry(param_pts);   % their images on boundary of domain
circle_pts = bcfun.values;       % corresponding pts on unit circle 
npts_inner = 30;                 % map interior points too - TO BE CHANGED
ni2 = 2*npts_inner;
cc = 2*rand(ni2,1)+2i*rand(ni2,1)-(1+1i);
cc = cc(abs(cc)<1);
cc(1) = 0;
circle_inner = cc(1:npts_inner);
dd = zeros(size(circle_inner));
dsdt = abs(diff(bndry));
s = cumsum(dsdt); scl = 2*pi/s(2*pi);
s = scl*s; dsdt = scl*dsdt;      % arc length scaled to [0 2pi]
dzetads = diff(bcfun);
dzetads = dzetads(s);
bcfuns = bcfun(s);
bdd = bndry*dzetads*dsdt;
for k = 1:length(circle_inner);
   dd(k) = (1/(2i*pi))*sum(bdd/(bcfuns-circle_inner(k)));
end
domain_inner = dd;
%si = inv(s); domain_pts = bndry(si(param_pts));
domain_pts = [domain_pts; domain_inner]; domain_pts = domain_inner;
circle_pts = [circle_pts; circle_inner]; circle_pts = circle_inner;

%% Construct function handles for the map and its inverse
tol = 1e-10;
[f,pol] = aaa(circle_pts, domain_pts, 'tol', tol);
[finv,poli] = aaa(domain_pts, circle_pts, 'tol', tol);
if min(abs(poli))<1, disp('Warning! -- poles in disk'), end

%% Plot results
if nargin > 2
  err = norm(domain_pts-finv(f(domain_pts)),inf);
  t = chebfun('t', [0 2*pi]);
  FS = 'fontsize'; MS = 'markersize';
  HA = 'horizontalalignment'; RT = 'right';
  clf

  % Plot the problem domain:
  subplot(1,2,1)
  plot(bndry,'k'), hold on
  xx = chebpts(100); cc = exp(1i*pi*linspace(-1,1,300));
  for theta = (0:7)*pi/8
    z = exp(1i*theta)*xx; plot(finv(z),'k')
  end
  for r = (1:5)/6
    z = r*cc; plot(finv(z),'k')
  end
  %plot(pol,'.r',MS,12)
  plot(domain_pts,'.b',MS,6)
  axis square, axis(1.3*[-1 1 -1 1])
  set(gca,FS,8)
  s1 = sprintf('f: type (%d,%d)', length(pol), length(pol));
  s2 = sprintf('finv: type (%d,%d)', length(poli), length(poli));
  s3 = sprintf('max err: %4.1e',err);

  % Plot the disk:
  subplot(1,2,2)
  plot(exp(1i*t),'k'), hold on
  for theta = (0:7)*pi/8
    z = exp(1i*theta)*xx; plot(z+1e-8i,'k')
  end
  for r = (1:5)/6
    z = r*cc; plot(z,'k')
  end
  %plot(poli,'.r',MS,12)
  plot(circle_pts,'.b',MS,6)
  axis square, axis(1.3*[-1 1 -1 1])
  set(gca,FS,8)
  s1 = sprintf('f: type (%d,%d)', length(pol), length(pol));
  s2 = sprintf('finv: type (%d,%d)', length(poli), length(poli));
  s3 = sprintf('max err: %4.1e',err);
  title([s1 '   ' s2 '   ' s3],FS,10,HA,RT)

  disp('press Enter to see 10,000 dots'), pause
  ndots = 10000;
  zz = rand(2*ndots,1) + 1i*rand(2*ndots,1);
  zz = 2*zz-(1+1i);
  zz = zz(abs(zz)<1); zz = zz(1:ndots);
  tic, ww = finv(zz); zz2 = f(ww); t = toc;
  disp(['time to map them back and forth: ', num2str(t) ' secs.'])
  err = norm(zz2-zz,inf);
  disp(['error in these points: ', num2str(err)])
  plot(zz,'.k',MS,3)
  subplot(1,2,1), plot(ww,'.k',MS,3)

end
