function phi_A = phi_eval(A,npts)

% uses pieces of conformal.m and kerzstein.m to evaluate phi(A), where phi
% is the conformal map from W(A) -> D and npts is the number of points used
% along the boundary of W(A) in the cauchy integral discretization

W = fov(A);
L = fov(A');

n = length(A);
ctr = trace(A)/n;

bndry = chebfun(conj(L),[0 2*pi],'trig');

dbndry = diff(bndry);
s = cumsum(abs(dbndry));
S = s(end);

% dbndry_arg = angle(dbndry);
% dbndry_arg = unwrap(dbndry_arg-dbndry_arg(0));
% clockwise = (diff(dbndry_arg([0 end])) < 0);
% % This is where the time gets spent! :
% if clockwise
%   u = inv(S-s); 
% else
%   u = inv(s);
% end
u = inv(s);

minu = u(0); maxu = u(end);
bndry2 = newDomain(bndry,[minu maxu]);
W = bndry2(u);
Wprime = diff(W);

ds = S/npts;
svec = [0:npts-1]'*ds;
Wvec = W(svec);

z = Wvec;
dzds = Wprime(svec);
phi_z = kerzstein(bndry, W, Wprime, S, npts, ctr); % compute images on D

phi_A = zeros(n);

for j=1:npts
    phi_A = phi_A + phi_z(j)*dzds(j)*ds*eye(n)/(z(j)*eye(n) - A);
end

phi_A = phi_A/(2i*pi);

end