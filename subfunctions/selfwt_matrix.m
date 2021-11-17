function selfwt=selfwt_matrix(elemType,normal_order,gamma,node,element)


% % Generates the force matrix due to self weight
% 
% numelem=size(element,1);
% numnode = size(node,1);
% total_unknown = numnode*2;
% selfwt=zeros(total_unknown,1);
% 
% for iel = 1 : numelem 
%     sctr1 = element(iel,:);    % element connectivity    
%     swpt=sctr1.*2;             %element degree of freedom
%    
%     [W,Q] = gauss_rule(iel,elemType,normal_order);                      
%     
%     for q=1:size(W,1)
%         pt = Q(q,:);
%         wt = W(q);                           % quadrature point
%         [N,dNdxi]=shape_func(elemType,pt);
%         J0 = node(sctr1,:)'*dNdxi;          % element Jacobian matrix
%         
%         selfwt(swpt) = selfwt(swpt)+ N*(-1*gamma)*det(J0)*wt;
%     end
%     
% end    
% 
% end   % end of function

% Generates the force matrix due to self weight

numelem=size(element,1);
numnode = size(node,1);
total_unknown = numnode*2;
selfwt=zeros(total_unknown,1);
k1=1;
for iel = 1 : numelem 
    sctr1 = element(iel,:);    % element connectivity    
    swpt=sctr1.*2;             %element degree of freedom
   
    [W,Q] = gauss_rule(iel,elemType,normal_order);                      
    
    for q=1:size(W,1)
        pt = Q(q,:);
        wt = W(q);                           % quadrature point
        [N,dNdxi]=shape_func(elemType,pt);
        J0 = node(sctr1,:)'*dNdxi;          % element Jacobian matrix
        %%%%% Additional/modification
          y_cord  = node(sctr1(:,2),2);
          r1 = size(gamma,2);
               for z1 = 1:r1
                   if y_cord >= gamma(z1,1) && y_cord < gamma(z1,2)
                    gamma_eq = gamma(z1,3);   
                    break;
                   end
               end
        selfwt(swpt) = selfwt(swpt)+ N*(-1*gamma_eq)*det(J0)*wt;
        temp1(k1)= gamma_eq; %% for varifications
        k1=k1+1;
    end
    
end    
varify1 = unique(temp1);
end   % end of function

