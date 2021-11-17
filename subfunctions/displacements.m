function [U,u_x,u_y] =displacements(dispNodes,dispNodes1,numnode,K,f,selfwt)

% evaluates the unknown degree of freedom (displacements) at the nodes

total_unknown=2*numnode;

udofs  = [(dispNodes.*2)-1;(dispNodes1.*2)-1]; %prescribed disp.in x-dir
vdofs  = dispNodes.*2;                         %prescribed disp. in y-dir
dofs=union(udofs(:),vdofs(:));                 %overall prescribed disp.
unknowndof=setdiff((1:total_unknown)',dofs);

F=f(unknowndof)+selfwt(unknowndof);
u=K(unknowndof,unknowndof)\F;
U=zeros(total_unknown,1);
U(unknowndof)=u;
u_x = U(1:2:2*numnode-1) ; % 1 3 5 7 ...
u_y = U(2:2:2*numnode) ; % 2 4 6 8 ...

end    % end of function

