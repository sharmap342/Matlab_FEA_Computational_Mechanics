function [node,element] =mesh_t6_elem(L,D,numx,numy)
%Forms the element and node matrices for T6 element, with the element
%and node matrices arranged in a nodal counterclockwise order 
%(primary nodes then secondary nodes).

 xx=repmat((0:L/(2*numx):L)',2*numy+1,1);
 yy=sort(repmat((-1*D/2:D/(2*numy):D/2)',2*numx+1,1));
 node=[xx yy];
 inc_u=[2 2 2 2 2 2];
 inc_v=[2+(4*numx) 2+(4*numx) 2+(4*numx) 2+(4*numx) 2+(4*numx) 2+(4*numx)];
     node_pattern1=[1 2 3 2*numx+3 4*numx+3 2*numx+2];
     node_pattern2=[3 2*numx+4 4*numx+5 4*numx+4 4*numx+3 2*numx+3];
     elements1=make_elem(node_pattern1,numx,numy,inc_u,inc_v);
     elements2=make_elem(node_pattern2,numx,numy,inc_u,inc_v);
           
      numberelem=2*numx*numy;
      elements=zeros(numberelem,6);
      elements(1:2:end,:)=elements1;
      elements(2:2:end,:)=elements2;
          element(:,1)=elements(:,1);
          element(:,2)=elements(:,3);
          element(:,3)=elements(:,5);
          element(:,4)=elements(:,2);
          element(:,5)=elements(:,4);
          element(:,6)=elements(:,6);
end

