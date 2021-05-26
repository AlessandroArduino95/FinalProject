%% PARABOLIC PDEs SOLVER FOR THE COURSE "NUMERICAL METHODS"
%  BY Alessandro Arduino

%% CLEAN UP

clear all
close all
clc

%% DATA INPUT

min = -1;
max = 1;
step = 1;

%% Setup of nodes and mesh

x = min:step:max; % calculate x coordinates
y = min:step:max; % calculate y coordinates 

edge = length(x); % meaning the number of nodes on one side

[x,y] = meshgrid(x,y); % combinations of x and y

x = x(:); % convert matrix of x into vector of x coordinates
y = y(:); % convert matrix of y into vector of y coordinates


nodes = [x y]; % node matrix (each line is a node)
numnodes = edge^2; % total number of nodes

H = zeros(numnodes,numnodes);  % preallocation of the H matrix (size = number of nodes)
P = zeros(numnodes,numnodes);  % preallocation of the P matrix (size = number of nodes)
N1 = zeros(3,3); % preallocation of the N1 matrix for the nodes for the calculations needed in the next step
N2 = zeros(3,3); % preallocation of the N2 matrix for the nodes for the calculations needed in the next step
N3 = zeros(3,3); % preallocation of the N3 matrix for the nodes for the calculations needed in the next step

% The mesh is made so that the elements are numbered from bottom left to
% the top right as are the nodes suche as in the below model
%
%  3-------4
%  |  1   /|
%  |    /  |
%  |  /    |
%  |/   2  |
%  1-------2
%
% In the following, the nodes will be taken in order in order to calculate
% the hij of the elements 
%% Calculate H nodes

for count = 1:numnodes % did I check all the nodes?
    
    if (mod(count,edge) ~= 0 && count <= edge*(edge-1)) % skipping every node at the end of the line and the top line (the equal is actually not necessary as it skipped by the inequality)
    
    % FIRST ELEMENT ASSOCIATED WITH NODE count
    %
    %   3-------2
    %   |  1   /
    %   |    /  
    %   |  /    
    %   |/     
    %   1
    %
    node1 = [1,nodes(count,:)]; %I'm looking at the two elements with node1 as the lower left corner
    node2 = node1 + [0,1,1]; % is diagonal to the first
    node3 = node2 + [0,-1,0]; % above the first
    
    % coordinates needed for the P matrix
    x1 = node1(1,2);  
    x2 = node2(1,2);
    y2 = node2(1,3);
    
    % first element with node node 1 using each point and moving
    % counterclockwise (assures the same determinant for all matrices)
    N1 = [node1; node2; node3]; 
    N2 = [node2; node3; node1];
    N3 = [node3; node1; node2];
    
    index1 = count; % I'm looking at node #count
    index3 = count + edge; % The node 3 is right above node 1
    index2 = index3 +1; % Node 2 is 1 step to the right of node 2
    
    % Calculate the parameters related to each 
    
    p1 = N1\[1;0;0];
    p2 = N2\[1;0;0];
    p3 = N3\[1;0;0]; 
    
    A = det(N1)/2; % as the mesh is regular, the determinant is the same everywhere, it would be possible to calculate it outside the loop but the imrpovement in speed wouldn't be significant yhus is not done
    
    % Calculate the hij values for the three nodes and their combinations 
    
    H(index1,index1) = H(index1,index1) + (p1(2)*p1(2)+p1(3)*p1(3))/(4*A);
    H(index1,index2) = H(index1,index2) + (p1(2)*p2(2)+p1(3)*p2(3))/(4*A);
    H(index1,index3) = H(index1,index3) + (p1(2)*p3(2)+p1(3)*p3(3))/(4*A);
    H(index2,index2) = H(index2,index2) + (p2(2)*p2(2)+p2(3)*p2(3))/(4*A);
    H(index2,index3) = H(index2,index3) + (p2(2)*p3(2)+p2(3)*p3(3))/(4*A);
    H(index3,index3) = H(index3,index3) + (p3(2)*p3(2)+p3(3)*p3(3))/(4*A);
    
    % Update the the values on the other side of diagonal (ensures
    % symmetry)
    
    H(index2,index1) = H(index1,index2);
    H(index3,index1) = H(index1,index3);
    H(index3,index2) = H(index2,index3);
    
        % Calculate the pij values
    P(index1,index1) = P(index1,index1) + pij_upp(x1,x2,y2,p1,p1);
    P(index1,index2) = P(index1,index2) + pij_upp(x1,x2,y2,p1,p2);
    P(index1,index3) = P(index1,index3) + pij_upp(x1,x2,y2,p1,p3);
    P(index2,index2) = P(index2,index2) + pij_upp(x1,x2,y2,p2,p2);
    P(index2,index3) = P(index2,index3) + pij_upp(x1,x2,y2,p2,p3);
    P(index3,index3) = P(index3,index3) + pij_upp(x1,x2,y2,p3,p3);
    
    % Update the the values on the other side of diagonal (ensures
    % symmetry)
    
    P(index2,index1) = P(index1,index2);
    P(index3,index1) = P(index1,index3);
    P(index3,index2) = P(index2,index3);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%------------------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%----------------------%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SECOND ELEMENT ASSOCIATED WITH NODE count
    %
    %           3
    %          /|
    %        /  |
    %      /    |
    %    /      |
    %   1-------2
    %
    node1 = [1,nodes(count,:)]; %I'm looking at the two elements with node1 as the lower left corner
    node2 = node1 + [0,1,0]; % to the right of the first one
    node3 = node2 + [0,0,1]; % above the second one
    
    % coordinates needed for the P matrix
    x1 = node1(1,2);  
    x2 = node2(1,2);
    y1 = node1(1,3);
    
    % second element with node node 1 using each point and moving
    % counterclockwise (assures the same determinant for all matrices)
    N1 = [node1; node2; node3]; 
    N2 = [node2; node3; node1];
    N3 = [node3; node1; node2];
    
    index1 = count; % I'm looking at node #count
    index2 = count + 1; % The node 2 is to the right of node 1
    index3 = index2 + edge; % Node 3 atop node 2                 
    
    % Calculate the parameters related to each 
    
    p1 = N1\[1;0;0];
    p2 = N2\[1;0;0];
    p3 = N3\[1;0;0]; 
        
    % Calculate the hij values for the three nodes and their combinations 
    
    H(index1,index1) = H(index1,index1) + (p1(2)*p1(2)+p1(3)*p1(3))/(4*A);
    H(index1,index2) = H(index1,index2) + (p1(2)*p2(2)+p1(3)*p2(3))/(4*A);
    H(index1,index3) = H(index1,index3) + (p1(2)*p3(2)+p1(3)*p3(3))/(4*A);
    H(index2,index2) = H(index2,index2) + (p2(2)*p2(2)+p2(3)*p2(3))/(4*A);
    H(index2,index3) = H(index2,index3) + (p2(2)*p3(2)+p2(3)*p3(3))/(4*A);
    H(index3,index3) = H(index3,index3) + (p3(2)*p3(2)+p3(3)*p3(3))/(4*A);
    
    % Update the the values on the other side of diagonal (ensures
    % symmetry)
    
    H(index2,index1) = H(index1,index2);
    H(index3,index1) = H(index1,index3);
    H(index3,index2) = H(index2,index3);
    
    % Calculate the pij values
    P(index1,index1) = P(index1,index1) + pij_low(x1,x2,y1,p1,p1);
    P(index1,index2) = P(index1,index2) + pij_low(x1,x2,y1,p1,p2);
    P(index1,index3) = P(index1,index3) + pij_low(x1,x2,y1,p1,p3);
    P(index2,index2) = P(index2,index2) + pij_low(x1,x2,y1,p2,p2);
    P(index2,index3) = P(index2,index3) + pij_low(x1,x2,y1,p2,p3);
    P(index3,index3) = P(index3,index3) + pij_low(x1,x2,y1,p3,p3);
    
    % Update the the values on the other side of diagonal (ensures
    % symmetry)
    
    P(index2,index1) = P(index1,index2);
    P(index3,index1) = P(index1,index3);
    P(index3,index2) = P(index2,index3);
    end
end
