function [usol,rf,kf,rk] =PCG_Jacobi(H,f,epsilon,kmax,u0)
%PCG function to solve the linear system

%Input
% Linear system Hu=f
% relative tolerance epsilon
% Maximum iteration number kmax

% Output
% Approximate solution usol
% Final relative residual norm: rf
% Number of iteration: kf
% Vector of relative residual norm for each iteration: rk


M_inv= sparse(diag(1./diag(H))); % Jacobi preconditioner M^-1

%CGM initialization
u=u0; 
r=f-H*u; 
tau=norm(r); %Norm of the k-th residual 
p=M_inv*r ; %Research direction

%CGM iteration
k=0; 
Rho=r'*M_inv*r;
tauk = 0;

while tau>epsilon*norm(f) && k<kmax
z=H*p;
alpha=p'*r/(p'*z);
u=u+alpha*p;
r=r-alpha*z;
g=M_inv*r;
Rho_p=Rho;
Rho=r'*g;
beta =Rho/Rho_p;
p=g+beta*p;
tau=norm(r);
tauk(k+1)=tau;
k=k+1;
end
% Solution
usol=u;
rf=tau/norm(f);
kf=k;
rk=tauk/norm(f);
end