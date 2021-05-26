function [f,H] = dirichletBoundary(H,f,xcoord,ycoord,u)
    % this codes automatically sets the dirichlet boundary conditions at
    % the domain's upper and lower edges.
    
    edgex = length(xcoord);
    edgey = length(ycoord);
    upleft = edgex*(edgey-1) + 1;
    dirNodes = [1:1:edgex, upleft:1:(upleft+edgex-1) ];
    penalty = 10^20;


    for i = 1:length(dirNodes)
        node = dirNodes(i);
        H(node,node) = penalty;
        f(node) = u(node)*penalty;
    end
end