function [norm_B] = bp_norm(alpha,A)
% blaschke optimization objective function

% computes norm of B(z), where B is a degree n-1 Blaschke product of A 
% with complex roots given by alpha (note real and imaginary parts are
% split so that alpha is a 2*(n-1) length vector here

n = length(A);
B = eye(n);

r = alpha(1:n-1);
theta = alpha(n:end);

for j=1:n-1
    numer = (A-r(j)*exp(1i*theta(j))*eye(n));
    denom = (eye(n) - r(j)*exp(-1i*theta(j))*A);
    B = B*numer/denom;
end
norm_B = norm(B);

end
