function Bfem4 =Bmatrix4(pt,elemType,iel)
% Gives the strain displacement matrix (B matrix of size 4x8)of each element
global node element

sctr = element(iel,:);
nn   = length(sctr);
[N,dNdxi] = shape_func(elemType,pt); % element shape functions
J0 = node(sctr,:)'*dNdxi;             % element Jacobian matrix
invJ0 = inv(J0);
dNdx  = dNdxi*invJ0;                  % derivatives of N w.r.t XY
Gpt = N'*node(sctr,1);              % GP in global coord, used

    Bfem4 = zeros(4,2*nn);
    Bfem4(1,1:2:2*nn)  = dNdx(:,1)' ;
    Bfem4(2,2:2:2*nn)  = dNdx(:,2)' ;
    Bfem4(3,1:2:2*nn)  = dNdx(:,2)' ;
    Bfem4(3,2:2:2*nn)  = dNdx(:,1)' ;
    

end  % end of function                    


    