function [u_next] = thetaMethod(H,P,u_prev,f,theta,dt,epsilon,kmax,xcoord,ycoord,ustationary,q_fun,nodes)
    % CODE TO CALCULATE THE NEXT TIME SOLUTION
    K1 = theta*H + P*(1/dt);
    K2 = P*(1/dt) -(1-theta)*H;
    
    ftheta = K2*u_prev -theta*f -(1-theta)*f;
    
    u0 = ones(length(u_prev),1);
    
    [ftheta,K1] = dirichletBoundary(K1,ftheta,xcoord,ycoord,ustationary); % I can use the solution at steady state as I'm only interested in the fixed values on the boundary
    ftheta = neumanntBoundary(ftheta,xcoord,ycoord,q_fun,nodes);
    ftheta=ftheta';
    [usol,~,~,~] = PCG_Jacobi(K1,ftheta,epsilon,kmax,u0); % solve the linear system


    u_next = usol;


end