% optimal blaschke product calculation

% this script does the following:

% 1. defines a matrix A whose optimal blaschke product we wish to find
% 2. use conformal code to map phi: W(A) -> unit disk
% 3. evaluate phi(A)
%    a. diagonalizable and well-conditioned
%    b. non-diagonalizable or poorly conditioned (cauchy integral)
% 4. use fmincon to find optimal blaschke product for A (i.e. phi(A))
% 5. the optimal blaschke product is returned as B, with roots alpha

%% compute conformal map phi:W(A) -> D, form phi(A)

% input matrix A here

B = [-2 7; 0 -2];
A = (B+eye(2))/(B - eye(2));

W = fov(A);
L = fov(A');

n = length(A);
ctr = trace(A)/n;

% ctr = 0;
% ctr = mean(W);

bndry = chebfun(conj(L),[0 2*pi],'trig');
[f, finv] = conformal(bndry, ctr);

% assuming A is diagonalizable 

% [V,D] = eig(A);
% phi_A = V*diag(f(diag(D)))/V;

% if A is not diagonalizable, need to compute via cauchy integral formula

phi_A = phi_eval(A,2^8); 

%% solve for blaschke product via fmincon

norm_phi_A = norm(phi_A);
a0 = zeros(1,2*n-2);

LB = zeros(1,2*n-2);
UB = ones(1,2*n-2);
UB(n:2*n-2) = 2*pi;

[a, B_norm] = fmincon(@(a) -bp_norm(a,phi_A),a0,[],[],[],[],LB,UB);
B_norm = -B_norm;

optimal_a = a;
optimal_B_norm = B_norm;

n_tries = 10;

for j = 1:n_tries
    a0 = rand(1,2*n-2);
    a0(n:2*n-2) = a0(n:2*n-2)*2*pi;
    
    [a, B_norm] = fmincon(@(a) -bp_norm(a,phi_A),a0,[],[],[],[],LB,UB);
    B_norm = -B_norm;
    
    if B_norm > optimal_B_norm
        optimal_a = a;
        optimal_B_norm = B_norm;
    end   
end

alpha = optimal_a(1:n-1).*exp(1i*optimal_a(n:2*n-2));

rr = optimal_a(1:n-1);
theta = optimal_a(n:end);

B = eye(n);

for j=1:n-1
    numer = (phi_A-rr(j)*exp(1i*theta(j))*eye(n));
    denom = (eye(n) - rr(j)*exp(-1i*theta(j))*phi_A);
    B = B*numer/denom;
end

optimal_B_norm;
alpha;
abs(alpha);

if optimal_B_norm < norm(phi_A)
    optimal_B_norm = norm(phi_A);
    alpha = [1 0];
end

if optimal_B_norm < norm(phi_A^2)
    optimal_B_norm = norm(phi_A^2);
    alpha = [0 0];
end

optimal_B_norm  % norm of optimal blaschke product
alpha           % root of optimal blaschke product
abs(alpha)      % norms of roots
%plot(W)
