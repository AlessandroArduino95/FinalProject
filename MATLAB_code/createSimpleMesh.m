function [coordinates,topology,S_vect] = createSimpleMesh(x_val,y_val,S)
% Creates the matrix containing the coordinates and topology of a simple,
% regular mesh of a rectangular region


[y_vect,x_vect] = meshgrid(x_val,y_val); % combinations of x and y


x_vect = x_vect(:); % convert matrix of x into vector of x coordinates
y_vect = y_vect(:); % convert matrix of y into vector of y coordinates

edgex = length(x_val); % number of nodes on the x side
edgey = length(y_val); % number of nodes on the y side

numel = (edgex - 1)*(edgey - 1)*2;

coordinates = [x_vect y_vect]; % node matrix (each line is a node)
numnodes = size(coordinates,1);


temp_top = zeros(numel,3);
i = 1;
for node = 1:numnodes % did I check all the nodes?

            if (mod(node,edgex) ~= 0 && node <= edgex*(edgey-1)) % skipping every node at the end of the line and the top line (the equal is actually not necessary as it skipped by the inequality)
                
                temp_top(i,:) = [ node, node + edgex + 1, node + edgex];
                
                temp_top(i+1,:) = [ node, node + 1, node + 1 + edgex];
                
                i = i + 2;
            end
end

topology = temp_top;
S_vect = S*ones(length(topology)); % creates a vector containing the values of S for all the element, in this case it is fixed but it could actually be loaded according to the real situation

end