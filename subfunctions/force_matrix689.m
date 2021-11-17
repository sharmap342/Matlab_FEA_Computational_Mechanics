function [f,sctry]=force_matrix689(node,topEdge,topEdge1,sigmatoy,sigmatox,load_edge1,load_edge2)

% Generates the force matrix due to externally applied loads 
% for Q9,Q8 and T6 elements and the location of these loads
numnode = size(node,1);
total_unknown = numnode*2;
f = zeros(total_unknown,1);

[W,Q]=gauss_pt_wt(2,'GAUSS',1);

ii1=intersect(find(node(:,1)==load_edge1),unique(topEdge));
jj1=intersect(find(node(:,1)==load_edge2),unique(topEdge));
ii2=find(topEdge(:,1)==ii1);
jj2=find(topEdge(:,2)==jj1);
ee2=topEdge1(1:end,2).*2;ee4=ee2-1;
ee1=unique(topEdge); 
ee3=ee1(1:2:end).*2;ee6=ee3-1;

for e =ii2:jj2
        sctr = topEdge(e,:);
        sctry = sctr.*2 ;
        sctrx = sctr.*2-1;
    for q=1:size(W,1)
        pt = Q(q,:);
        wt = W(q);
        [N,dNdxi]  = shape_func('L2',pt);
        J0 = abs( node(sctr(2))-node(sctr(1)) )/2;
        f(sctry) = f(sctry) + N*sigmatoy*det(J0)*wt;
        f(sctrx) = f(sctrx) + N*sigmatox*det(J0)*wt;
    end   % of quadrature loop
end       % of element loop
f(ee3)=f(ee3)*(2/3);
f(ee2)=f(ee2)*(4/3);
f(ee6)=f(ee6)*(2/3);
f(ee4)=f(ee4)*(4/3);
end   % end of function
