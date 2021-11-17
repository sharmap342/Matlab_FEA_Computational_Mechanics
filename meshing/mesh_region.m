function [node,element] = mesh_region(pt1, pt2, pt3, pte1, pt4,numx,numy,elemType)

    % Generates an array of nodal connectivity (coordinates of each node)
    % and element connectivity (nodes of each element with counterclockwise
    % ordering.)
    % given the four corners points of the domain,number of elements
    % in each direction (numx,numy)and the element type (Q4,Q8,Q9 and T3)
    
global L D
nnx = numx+1;
nny = numy+1;

switch elemType
        
    case 'Q4'           
        node=square_node_array(pt1,pt2,pt3,pte1,pt4,nnx,nny);
        inc_u=1;
        inc_v=nnx;
        node_pattern=[ 1 2 nnx+2 nnx+1 ];
        [element]=make_elem(node_pattern,numx,numy,inc_u,inc_v);

    case 'Q9'           
        [node,element]=structured_q9_mesh(pt1,pt2,pt3,pt4,numx,numy);
        
    case 'Q8'
        [node,element]=structured_q8_mesh(pt1,pt2,pt3,pt4,numx,numy);
        
    case 'T3'           
        
        node=square_node_array(pt1,pt2,pt3,pt4,nnx,nny);
        node_pattern1=[ 1 2 nnx+1 ];
        node_pattern2=[ 2 nnx+2 nnx+1 ];
        inc_u=1;
        inc_v=nnx;
        numberelem=2*numx*numy;
        element=zeros(numberelem,3);
        element1=make_elem(node_pattern1,numx,numy,inc_u,inc_v);
         element2=make_elem(node_pattern2,numx,numy,inc_u,inc_v);
       element(1:2:end,:)=element1;
       element(2:2:end,:)=element2;
    otherwise
        error('only Q4,Q9,Q8,and T3 are supported by this function');
end
end



