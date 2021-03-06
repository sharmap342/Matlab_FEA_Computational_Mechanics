function eldof =elementdof(elemType,sctr)

% forms the global degree of freedom nodes in each element from the node
% identification number.

switch (elemType)
       case 'Q9'
        eldof=[sctr(1)*2-1;sctr(1)*2;...
               sctr(2)*2-1;sctr(2)*2;...
               sctr(3)*2-1;sctr(3)*2;...
               sctr(4)*2-1;sctr(4)*2;...
               sctr(5)*2-1;sctr(5)*2;...
               sctr(6)*2-1;sctr(6)*2;...
               sctr(7)*2-1;sctr(7)*2;...
               sctr(8)*2-1;sctr(8)*2;...
               sctr(9)*2-1;sctr(9)*2];
       case 'Q8'
        eldof=[sctr(1)*2-1;sctr(1)*2;...
               sctr(2)*2-1;sctr(2)*2;...
               sctr(3)*2-1;sctr(3)*2;...
               sctr(4)*2-1;sctr(4)*2;...
               sctr(5)*2-1;sctr(5)*2;...
               sctr(6)*2-1;sctr(6)*2;...
               sctr(7)*2-1;sctr(7)*2;...
               sctr(8)*2-1;sctr(8)*2];
       case 'Q4'
        eldof=[sctr(1)*2-1;sctr(1)*2;...
               sctr(2)*2-1;sctr(2)*2;...
               sctr(3)*2-1;sctr(3)*2;...
               sctr(4)*2-1;sctr(4)*2]; 
       case 'T6'
        eldof=[sctr(1)*2-1;sctr(1)*2;...
               sctr(2)*2-1;sctr(2)*2;...
               sctr(3)*2-1;sctr(3)*2;...
               sctr(4)*2-1;sctr(4)*2;...
               sctr(5)*2-1;sctr(5)*2;...
               sctr(6)*2-1;sctr(6)*2];
       case 'T3'
         eldof=[sctr(1)*2-1;sctr(1)*2;...
                sctr(2)*2-1;sctr(2)*2;...
                sctr(3)*2-1;sctr(3)*2]; 
       case 'L3'
         eldof=[sctr(1)*2-1;sctr(1)*2;...
               sctr(2)*2-1;sctr(2)*2;...
               sctr(3)*2-1;sctr(3)*2]; 
       case 'L2'
         eldof=[sctr(1)*2-1;sctr(1)*2;...
               sctr(2)*2-1;sctr(2)*2];
               
end

end   % end of function

