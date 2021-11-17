function Bfem =Bmatrix(pt,elemType,iel)
% Gives the strain displacement matrix (B matrix of size 3x8)of each element
global node element

sctr = element(iel,:);
nn   = length(sctr);
[N,dNdxi] = shape_func(elemType,pt);  % element shape functions
J0 = node(sctr,:)'*dNdxi;             % element Jacobian matrix
invJ0 = inv(J0);
dNdx  = dNdxi*invJ0;                  % derivatives of N w.r.t XY
Gpt = N'*node(sctr,:);               % GP in global coord, used


    Bfem = zeros(3,2*nn);
    Bfem(1,1:2:2*nn)  = dNdx(:,1)' ;
    Bfem(2,2:2:2*nn)  = dNdx(:,2)' ;
    Bfem(3,1:2:2*nn)  = dNdx(:,2)' ;
    Bfem(3,2:2:2*nn)  = dNdx(:,1)' ;

end  % end of function                    


    