function u_next = forwardDifferenceMethod(H,P,u_prev,f,dt,epsilon,kmax)
    % CODE TO CALCULATE THE NEXT TIME SOLUTION
    K1 = P*(1/dt);
    K2 = P*(1/dt) - H;
    
    ftheta = K2*u_prev + f;
    
    u0 = zeros(length(u_prev),1);
    
    [usol,~,~,~] = PCG_Jacobi(K1,ftheta,epsilon,kmax,u0);


    u_next = usol;


end

