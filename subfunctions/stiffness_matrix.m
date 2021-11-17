function K=stiffness_matrix(node,element,elemType,normal_order,C)

%  Generates the element and global stiffness matrix

numnode = size(node,1);
numelem = size(element,1);
total_unknown = numnode*2;            
K = zeros(total_unknown,total_unknown);            

for iel = 1 : numelem 
    sctr = element(iel,:);               % element connectivity
    nn   = length(sctr);                 % number of nodes per element    
    eldof =elementdof(elemType,sctr);    %element degree of freedom
    
    [W,Q] = gauss_rule(iel,elemType,normal_order);              
      
        
    for kk = 1 : size(W,1)             % Loop on Gauss points 
        pt = Q(kk,:);                  % quadrature point
       
        [N,dNdxi] = shape_func(elemType,pt);  
        J0 = node(sctr,:)'*dNdxi;                % element Jacobian matrix
        Bfem =Bmatrix(pt,elemType,iel);
       
       
       K(eldof,eldof) = K(eldof,eldof)+Bfem'*C(1:3,1:3)*Bfem*W(kk)*det(J0);
    end                  % end of looping on Gauss Points
end

end    % end of function

