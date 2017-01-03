function [g_B, pts, images, W, Wprime, S] = kerzstein(L, npts, a)
%
%  Given a smooth (chebfun) curve L, defining the boundary of a region Omega,
%  this routine finds the boundary correspondence function g_B for a conformal
%  mapping from Omega to the unit disk, mapping a to 0.  It first finds npts
%  points (pts) on L that are equally spaced in arclength.  It then solves
%  an integral equation to determine the images of these points (images) on
%  the unit circle.  Finally, a trigonometric chebfun g_B is defined to map
%  the arclengths associated with pts to the images on the unit circle.
%  It also returns W, which is the curve L expressed as a function
%  of arclength, and Wprime, which is the derivative of W with respect to 
%  arclength.  Additionally, it returns the total arclength S of the curve L.

%  Express L as a function of arclength and compute its derivative.

S = arcLength(L);
s = cumsum(abs(diff(L)));

signed_area = sum(real(L).*diff(imag(L)));
direction = sign(signed_area);

if direction == 1     % Reverse the order of L if it is clockwise.
    u = inv(s);
else
    u = inv(S-s);
end

L = newDomain(L,minandmax(u));
W = L(u);
Wprime = diff(W);

%  Set up and solve the integral equation.

ds = S/npts;
svec = (0:npts-1)'*ds;
Wvec = W(svec);
pts = Wvec;
Wprimevec = Wprime(svec);
gammadotvec = Wprimevec ./ abs(Wprimevec);   % Unit tangents
gvec = (1/(2*pi*1i)) * conj(gammadotvec ./ (a-Wvec));   % RHS 

M = eye(npts,npts);
for j=1:npts,
  w = Wvec(j);
  for i=1:npts,
    z = Wvec(i);
    if i ~= j,
      M(i,j) = M(i,j) - (1/(2*pi*1i))*(conj(gammadotvec(i)/(z-w)) + gammadotvec(j)/(w-z))*ds;
    end;
  end;
end;
fvec = M\gvec;

%  Compute boundary correspondence function g_B.

Rprimevec = fvec.^2;
Rvec = (1/1i)*gammadotvec.*(Rprimevec./abs(Rprimevec));
images = Rvec;
g_B = chebfun(Rvec,[0,S],'trig');
