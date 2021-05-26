function value = pij_upp(x1,x2,y2,pi,pj)
    %  
    %   THIS CODE IS MEANT FOR LOWER TRIANGULAR ELEMENTS SUCH AS THE ONE
    %   BELOW
    %   Each time the code is invoked, it will take two nodes of the
    %   element (eventually coincidental ones)
    % 
    %   3-------2
    %   |  1   /
    %   |    /  
    %   |  /    
    %   |/     
    %   1
    %
    %
    % In order to make the calculations for the pij elements easier as it
    % would be possible to make more mistakes if this was directly part of
    % the main code
    ai = pi(1);  % again, to render the code more readable, each parameter is extracted and appropriately named
    bi = pi(2);  % in order to increase the readability
    ci = pi(3);
    aj = pj(1);
    bj = pj(2);
    cj = pj(3);
    
    % The big parameter has been divided acording to the calculations made
    % in my nodes in 6 smaller formulas.
    temp1 = ai*aj*((x2-x1)*y2-((x2^2-x1^2)/2));
    temp2 = ((y2*(x2^2-x1^2)/2)-((x2^3-x1^3)/3))*(ai*bj + aj*bi);
    temp3 = (((y2^2)/2)-((x2^3-x1^3)/6))*(ai*cj+aj*ci);
    temp4 = (((y2^2)*(x2^2-x1^2)/4)-((x2^4-x1^4)/8))*(ci*bj+cj*bi);
    temp5 = ((((x2^3-x1^3)*y2)/3)-((x2^4-x1^4)/4))*bi*bj;
    temp6 = ((((y2^3)*(x2-x1))/3)-((x2^4-x1^4)/12))*ci*cj;
    
    % sum each component to get the final value and then return
    value = temp1+temp2+temp3+temp4+temp5+temp6;
end