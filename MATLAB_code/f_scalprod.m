function f = f_scalprod(f_vect,nodes,topology)
    % Function meant to calculate the scalar product between f (vector of
    % known constants) and the basis of the function space where the
    % solution lies
    %
    ftemp = zeros(size(nodes,1),1);
    
        for i = 1:size(topology,1) % cycle through the elements
        
            el = topology(i,:); % extract the nodal indexes from the topology
        
            % Each node is located at index el(i) with i in [1-3]
            ind1 = el(1);
            ind2 = el(2);
            ind3 = el(3);
        
            N1 = [  1 nodes(ind1,:);
                    1 nodes(ind2,:);
                    1 nodes(ind3,:)];
    
            A = det(N1)/2; % as the mesh is regular, the determinant is the same everywhere, it would be possible to calculate it outside the loop but the imrpovement in speed wouldn't be significant yhus is not done
    
            ftemp(ind1) = ftemp(ind1) + f_vect(ind1)*(A/3);
            ftemp(ind2) = ftemp(ind2) + f_vect(ind2)*(A/3);
            ftemp(ind3) = ftemp(ind3) + f_vect(ind3)*(A/3);
        
        end
        
        f = ftemp;
end
   
 
            