function plot_field(node,connect,elem_type,field)
  
% Forms a color dependent finite element mesh for plotting outputs 
  
if ( size(field) == size(connect) )
  elementalField=1;
else
  elementalField=0;
end
% fill node if needed
if (size(node,2) < 3)
   for c=size(node,2)+1:3
      node(:,c)=[zeros(size(node,1),1)];
   end
end

holdState=ishold;
hold on

% plot elements
if ( strcmp(elem_type,'Q9') )      % Q9 element
  ord=[1,5,2,6,3,7,4,8,1];
elseif ( strcmp(elem_type,'T3') )  % T3 element
  ord=[1,2,3,1];
elseif ( strcmp(elem_type,'T6') )  % T6 element
  ord=[1,4,2,5,3,6,1];
elseif ( strcmp(elem_type,'Q4') )  % Q4 element
  ord=[1,2,3,4,1];
elseif ( strcmp(elem_type,'Q8') )  % Q8 element
   ord=[1,5,2,6,3,7,4,8,1];
elseif ( strcmp(elem_type,'L2') )  % L2 element
  ord=[1,2];   
end

for e=1:size(connect,1)
  
   xpt=node(connect(e,ord),1);
   ypt=node(connect(e,ord),2);      
   zpt=node(connect(e,ord),3);
   
   if ( elementalField )
     fpt=field(e,ord);
   else
     fpt=field(connect(e,ord));
   end
   
   fill3(xpt,ypt,zpt,fpt)
end

shading interp
axis equal
      
if ( ~holdState )
  hold off
end
end  % end of fuction
