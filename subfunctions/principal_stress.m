function [s1,s2,s3] = principal_stress(kk,iel,p,q,theta)

%calcualtes the principal stresses from the stress invariants

s1=p(:,kk,iel)+q(:,kk,iel)*(2/3)*sind(theta(:,kk,iel)-120);

s2=p(:,kk,iel)+q(:,kk,iel)*(2/3)*sind(theta(:,kk,iel));

s3=p(:,kk,iel)+q(:,kk,iel)*(2/3)*sind(theta(:,kk,iel)+120);


end % end of function