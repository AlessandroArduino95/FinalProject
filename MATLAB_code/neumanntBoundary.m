function f = neumanntBoundary(f,xcoord,ycoord,q_fun,nodes)
    % this codes automatically sets the Neumann boundary conditions at
    % the domain's left and right edges.
    % as the domain is regular, each node has two boundary sides connected
    % to it and thus the values of qi would be multiplied by 1
    edgex = length(xcoord);
    edgey = length(ycoord);
    
    neuNodesleft = [(edgex+1):edgex:(edgex*(edgey-2)+1)];
    neuNodesright = [2*edgex:edgex:edgex*(edgey-1)];
    
    for i = 1:length(neuNodesleft)
        node = neuNodesleft(i);
        xi = nodes(node,1);
        yi = nodes(node,2);
        f(node) = f(node) + q_fun(xi,yi);
    end
    
    for i = 1:length(neuNodesright)
        node = neuNodesright(i);
        xi = nodes(node,1);
        yi = nodes(node,2);
        f(node) = f(node) - q_fun(xi,yi);
    end
    
    f=f';
end