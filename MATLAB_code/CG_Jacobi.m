function [usol] =CG_Jacobi(H,f,epsilon,kmax,u0)
    %CG function to solve the linear system

    %Input
    % Linear system Hu=f
    % relative tolerance epsilon
    % Maximum iteration number kmax

    % Output
    % Approximate solution usol
    % Final relative residual norm: rf
    % Number of iteration: kf
    % Vector of relative residual norm for each iteration: rk


    %CGM initialization
    u=u0; 
    r_prev=f-H*u; % residual 
    p=r_prev ;    %Research direction
    k=0;          %CGM iteration

    while norm(r_prev)>epsilon*norm(f) && k<kmax
        z=H*p;              % calcolo z
        alpha=(p'*r_prev)/(z'*p);   % calcolo alpha
        u=u+alpha*p;        % calcolo nuovo guess
        r_next=r_prev-alpha*z;        % calcolo nuovo residuo
        beta = (r_next'*r_next)/(r_prev'*r_prev);    % calcolo beta
        p=r_next + beta*p;         % calcolo nuovo p
        r_prev = r_next;
        k=k+1;              % aggiorno il counter per le iterazioni
    end

    % Solution
    usol=u;
end