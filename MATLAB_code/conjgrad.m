function [x, error, niter, flag] = conjgrad(A, x, b, maxit, tol)
flag = 0; niter = 0; bnrm2 = norm( b );
if ( bnrm2 == 0.0 ), bnrm2 = 1.0; end
r = b - A*x; error = norm( r ) / bnrm2; P = r;
if ( error < tol ) return, end
for niter = 1:maxit
z = P \ r; rho = (r'*z);
if niter > 1
beta = rho / rho1; p = z + beta*p;
else
p = z;
end
q = A*p; alpha = rho / (p'*q );
x = x + alpha * p; r = r - alpha*q;
error = norm( r ) / bnrm2;
if ( error <= tol ), break, end
rho1 = rho;
end
if ( error > tol ) flag = 1; end