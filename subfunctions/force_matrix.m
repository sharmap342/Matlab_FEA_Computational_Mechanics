function [f,sctry]=force_matrix(node,topEdge,sigmatoy,sigmatox,load_edge1,load_edge2)

% Generates the force matrix due to externally applied loads
% for Q4 and T3 elements and the location of these loads

numnode = size(node,1);
total_unknown = numnode*2;
f = zeros(total_unknown,1);

[W,Q]=gauss_pt_wt(1,'GAUSS',1);

ii1=intersect(find(abs(node(:,1)-load_edge1)<=0.000001),unique(topEdge));
jj1=intersect(find(abs(node(:,1)-load_edge2)<=0.000001),unique(topEdge));
ii2=find(topEdge(:,1)==ii1);
jj2=find(topEdge(:,2)==jj1);

for e =ii2:jj2
       sctr = topEdge(e,:);
       sctry = sctr.*2 ;
       sctrx = sctr.*2-1;
    for q=1:size(W,1)
        pt = Q(q,:);
        wt = W(q);
        N  = shape_func('L2',pt);
        J0 = abs( node(sctr(2))-node(sctr(1)) )/2;
        f(sctry) = f(sctry) + N*sigmatoy*det(J0)*wt;
        f(sctrx) = f(sctrx) + N*sigmatox*det(J0)*wt;
    end   % end of quadrature loop
end       % end of element loop

end  % end of function

