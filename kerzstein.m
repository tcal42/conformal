function images = kerzstein(bndry, W, Wprime, S, npts, ctr)
%
%  Input: 
%
%  bndry = smooth periodic complex chebfun mapping [0, 2*pi]
%          defining the boundary of a region Omega
%  ctr = point to be mapped to 0 in the unit disk
%  npts = numbers of points on the boundary for discretization
%
%  Output:
%
%  images = images on the unit circle of npts points on boundary, 
%           equispaced by arclength s

%%
%  This subfunction was adapted from a code written by Trevor Caldwell and
%  Anne Greenbaum in 2016.

%  Reference:  Kerzman and Trummer, "Numerical conformal mapping via the Szego
%  kernel," Journal of Computational and Applied Mathematics 14 (1986), 111-123.

%  Express bndry as a function of arclength and compute its derivative.

%  Set up and solve the integral equation.

ds = S/npts;
svec = [0:npts-1]'*ds;
Wvec = W(svec);
pts = Wvec;
Wprimevec = Wprime(svec);
gammadotvec = Wprimevec ./ abs(Wprimevec);   % Unit tangents
gvec = (1/(2*pi*1i)) * conj(gammadotvec ./ (ctr-Wvec));

M = eye(npts);
for j = 1:npts
  gdj = gammadotvec(j);
  zwi = 1./(Wvec-Wvec(j)); zwi(j) = 0;
  M(:,j) = M(:,j) - (ds/(2i*pi))*(conj(zwi.*gammadotvec) - zwi*gdj);
end
fvec = M\gvec;

Rprimevec = fvec.^2;
Rvec = (1/1i)*gammadotvec.*(Rprimevec./abs(Rprimevec));
images = Rvec;
