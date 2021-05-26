function [P,H] = PHcalc(nodes,topology,S_vect)
    
    col = reshape(topology.',1,[])'; % column values
    rows = repelem(1:1:size(topology,1),3)'; % row values
    values = zeros(size(col,1),1); % intialization vector of zeros
    
    H = sparse(rows,col,values); % preallocation of the sparse H matrix
    H = H'*H;
    P = H;  % preallocation of the P matrix (size = number of nodes)
    
    for i = 1:size(topology,1) % cycle through the elements
        
        el = topology(i,:); % extract the nodal indexes from the topology
        
        % Each node is located at index el(i) with i in [1-3]
        ind1 = el(1);
        ind2 = el(2);
        ind3 = el(3);
        
        N1 = [  1 nodes(ind1,:);
                1 nodes(ind2,:);
                1 nodes(ind3,:)];
        
        % Calculate the parameters related to each 
    
        p1 = inv(N1)*[1;0;0].*det(N1);
        p2 = inv(N1)*[0;1;0].*det(N1);
        p3 = inv(N1)*[0;0;1].*det(N1); 
    
        A = det(N1)/2; % as the mesh is regular, the determinant is the same everywhere, it would be possible to calculate it outside the loop but the imrpovement in speed wouldn't be significant yhus is not done
    
        % Calculate the hij values for the three nodes and their combinations 
    
        H(ind1,ind1) = H(ind1,ind1) + (p1(2)*p1(2)+p1(3)*p1(3))/(4*A);
        H(ind1,ind2) = H(ind1,ind2) + (p1(2)*p2(2)+p1(3)*p2(3))/(4*A);
        H(ind1,ind3) = H(ind1,ind3) + (p1(2)*p3(2)+p1(3)*p3(3))/(4*A);
        H(ind2,ind2) = H(ind2,ind2) + (p2(2)*p2(2)+p2(3)*p2(3))/(4*A);
        H(ind2,ind3) = H(ind2,ind3) + (p2(2)*p3(2)+p2(3)*p3(3))/(4*A);
        H(ind3,ind3) = H(ind3,ind3) + (p3(2)*p3(2)+p3(3)*p3(3))/(4*A); 
        
        % Update the the values on the other side of diagonal (ensures
        % symmetry)
    
        H(ind2,ind1) = H(ind1,ind2);
        H(ind3,ind1) = H(ind1,ind3);
        H(ind3,ind2) = H(ind2,ind3);
        
        % P matrix entries
% %         P(ind1,ind1) = P(ind1,ind1) + A/6;
% %         P(ind1,ind2) = P(ind1,ind2) + A/12;
% %         P(ind1,ind3) = P(ind1,ind3) + A/12;
% %         P(ind2,ind2) = P(ind2,ind2) + A/6;
% %         P(ind2,ind3) = P(ind2,ind3) + A/12;
% %         P(ind3,ind3) = P(ind3,ind3) + A/6; 
        
         P(ind1,ind1) = P(ind1,ind1) + S_vect(i)*A/6;
         P(ind1,ind2) = P(ind1,ind2) + S_vect(i)*A/12;
         P(ind1,ind3) = P(ind1,ind3) + S_vect(i)*A/12;
         P(ind2,ind2) = P(ind2,ind2) + S_vect(i)*A/6;
         P(ind2,ind3) = P(ind2,ind3) + S_vect(i)*A/12;
         P(ind3,ind3) = P(ind3,ind3) + S_vect(i)*A/6; 

        % Update the the values on the other side of diagonal (ensures
        % symmetry)
    
        P(ind2,ind1) = P(ind1,ind2);
        P(ind3,ind1) = P(ind1,ind3);
        P(ind3,ind2) = P(ind2,ind3);
    end
end
