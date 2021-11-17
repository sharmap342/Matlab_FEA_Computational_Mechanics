function [Elements,Nodes]=q4totq9(element,node,numx,numy)

% forms the element and node matrices for nine node rectangular element 
% given the elemet and node matrices of a four node rectangular element 
% and the number of elements in x and y direction

  nnx=numx+1;
  num_u=numx;
  num_v=numy;
  inc_u=[2 2 2 2 2 2 2 2 2];
  if numx==3
    inc_v=[14 14 14 14 14 14 14 14 14];
  
  else
    inc_v=[14+4*(numx-3) 14+4*(numx-3) 14+4*(numx-3) 14+4*(numx-3)...
      14+4*(numx-3) 14+4*(numx-3) 14+4*(numx-3) 14+4*(numx-3) 14+4*(numx-3)];
    
  end
  
 node_pattern=[ 1 2 3 2*nnx+2 4*nnx+1 4*nnx 4*nnx-1 2*nnx 2*nnx+1];
inc=zeros(1,size(node_pattern,2));
e=1;
elements=zeros(num_u*num_v,size(node_pattern,2));

for row=1:num_v
   for col=1:num_u
      elements(e,:)=node_pattern+inc;
      inc=inc+inc_u;
      e=e+1;
   end
   inc=row*inc_v;
end
Elements=elements;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numNode=size(unique(Elements),1);
Nodes=zeros(numNode,2);
numElement=numx*numy;
for i=1:numElement
  
       Nodes(Elements(i,1),1)=node(element(i,1),1);
       Nodes(Elements(i,1),2)=node(element(i,1),2);
       
       Nodes(Elements(i,2),1)=(node(element(i,1),1)+node(element(i,2),1))/2;
       Nodes(Elements(i,2),2)= node(element(i,1),2);
       
       Nodes(Elements(i,3),1)=node(element(i,2),1);
       Nodes(Elements(i,3),2)=node(element(i,2),2);
       
       Nodes(Elements(i,4),1)=node(element(i,2),1);
       Nodes(Elements(i,4),2)=(node(element(i,2),2)+node(element(i,3),2))/2;
       
       Nodes(Elements(i,5),1)=node(element(i,3),1);
       Nodes(Elements(i,5),2)=node(element(i,3),2);
       
       Nodes(Elements(i,6),1)=(node(element(i,3),1)+node(element(i,4),1))/2;
       Nodes(Elements(i,6),2)= node(element(i,3),2);
       
       Nodes(Elements(i,7),1)=node(element(i,4),1);
       Nodes(Elements(i,7),2)=node(element(i,4),2);
       
       Nodes(Elements(i,8),1)=node(element(i,4),1);
       Nodes(Elements(i,8),2)=(node(element(i,4),2)+node(element(i,1),2))/2;
      
       Nodes(Elements(i,9),1)=(node(element(i,1),1)+node(element(i,2),1))/2;
       Nodes(Elements(i,9),2)=(node(element(i,4),2)+node(element(i,1),2))/2;
end
    Nodes=Nodes;
    Elements(:,1)=elements(:,1);
    Elements(:,2)=elements(:,3);
    Elements(:,3)=elements(:,5);
    Elements(:,4)=elements(:,7);
    Elements(:,5)=elements(:,2);
    Elements(:,6)=elements(:,4);
    Elements(:,7)=elements(:,6);
    Elements(:,8)=elements(:,8);
    Elements(:,9)=elements(:,9);
end
    
