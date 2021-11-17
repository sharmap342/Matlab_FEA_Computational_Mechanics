function element=make_elem(node_pattern,numx,numy,inc_u,inc_v)

% creates a connectivity list of primary nodes in Q4 and T3 element

if ( nargin < 5 )
   disp('Not enough parameters specified for make_elem function')
end

inc=[zeros(1,size(node_pattern,2))];
e=1;
element=zeros(numx*numy,size(node_pattern,2));

for row=1:numy
   for col=1:numx
      element(e,:)=node_pattern+inc;
      inc=inc+inc_u;
      e=e+1;
   end
     inc=row*inc_v;
end
end
