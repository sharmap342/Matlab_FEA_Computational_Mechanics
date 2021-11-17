function [dg1,dg2,dg3] =formdg(tsi,q,theta)

% calculates the partial derivatives of the plastic potential with respect 
% to p, J2 and J3.

   dg1=sind(tsi);
if sind(theta)>0.49   % close to theta=30 corner,smoothen the curve with
     sw=1;            % triaxial compression case s1 > s2 = s3
     dg2=(0.25/q)*(3-sw*sind(tsi));
     dg3=0;   
       
elseif -1*sind(theta)>049  % close to theta=-30 corner,smoothen the curve
     sw=-1;               % with triaxial extension case s1 = s2 > s3
     dg2=(0.25/q)*(3-sw*sind(tsi));
     dg3=0;
else                        % all other cases
    dg2=(sqrt(3)*cosd(theta)/(2*q))*(1+(tand(theta)*tand(3*theta))+...
        ((sind(tsi)/sqrt(3))*(tand(3*theta)-tand(theta))));
    dg3=1.5*((sqrt(3)*sind(theta))+(sind(tsi)*cosd(theta))/(q^2*cosd(3*theta)));

end
end % end of function