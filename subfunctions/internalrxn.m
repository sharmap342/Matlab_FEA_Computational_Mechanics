function [r] = internalrxn( node,element,elemType,normal_order,stress)

%Generates the vector of nodal force reactions due to internal stress. it
%can be used out side the iteration loop if the stresses at each gauss
%point is given.

numnode = size(node,1);
numelem = size(element,1);
total_unknown = numnode*2;  
r=zeros(total_unknown,1);

for iel = 1 : numelem 
    sctr = element(iel,:);                  % element connectivity
    nn   = length(sctr);                    % number of nodes per element    
    eldof =elementdof(elemType,sctr);       %element degree of freedom
   [W,Q] = gauss_rule(iel,elemType,normal_order);   % determine GP                       
       
    for kk = 1 :nn
        pt =Q(kk,:);                           
       [N,dNdxi] = shape_func(elemType,pt);
        J0 = node(sctr,:)'*dNdxi;                 
        Bfem4 =Bmatrix4(pt,elemType,iel);
       r(eldof) =r(eldof)+Bfem4'*stress(:,kk,iel)*W(kk)*det(J0);               
    end       % end of looping on Gauss Points
end             % end of looping on elements
end             % end of function