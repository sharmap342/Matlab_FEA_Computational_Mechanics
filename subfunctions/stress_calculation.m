function [stress,strain] =stress_calculation(node,element,elemType,U,normal_order,C)

% calculates the element strains and stresses at the nodes
% in x, y and xy directions. 

numelem=size(element,1);nonelm=size(element,2);stress=zeros(4,nonelm,numelem);strain=zeros(4,nonelm,numelem);
switch elemType
    case 'Q9'
      stresspoints=[-1 -1;1 -1;1 1;-1 1;0 -1;1 0;0 1;-1 0;0 0];
    case 'Q8'
      stresspoints=[-1 -1;1 -1;1 1;-1 1;0 -1;1 0;0 1;-1 0];
    case 'Q4'
      stresspoints=[-1 -1;1 -1;1 1;-1 1];
    case 'T3'
      stresspoints=[0 0;1 0;0 1];  
    otherwise
      stresspoints=[0 0;1 0;0 1;0.5 0;0.5 0.5;0 0.5];
end
for iel = 1 : numelem 
    sctr = element(iel,:); % element connectivity
    nn   = length(sctr);   % number of nodes per element
    
   eldof =elementdof(elemType,sctr);  
   [W,Q] = gauss_rule(iel,elemType,normal_order);   
        
    for kk = 1:nn
        pt = Q(kk,:);              % quadrature point      
        [N,dNdxi] = shape_func(elemType,pt);  % element shape functions
        J0 = node(sctr,:)'*dNdxi;             % element Jacobian matrix       
        Bfem4 =Bmatrix4(pt,elemType,iel);           
        strain(:,kk,iel)=Bfem4*U(eldof);
        strain(4,:,:)=0; 
        stress(:,kk,iel)=C*strain(:,kk,iel);  
      
    end                  % end of looping on GPs
end                      % end of looping on elements
end   % end of function

