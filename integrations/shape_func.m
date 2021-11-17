function [N,dNdxi]=shape_func(type,coord,dim)
  
% Gives the shape function and its derivatives with respect to x and y   
    
  if ( nargin == 2 )
    dim=1;
  end
  
  switch type
   case 'L2'  
    %%%%%%%%%%%%%%%%%%%%% L2 TWO NODE LINE ELEMENT %%%%%%%%%%%%%%%%%%%%% 
    %  
    %    1---------2
    %
    if size(coord,2) < 1
      disp('Error coordinate needed for the L2 element')
    else
     xi=coord(1);
     N=([1-xi,1+xi]/2)';
     dNdxi=[-1;1]/2;
    end  
    
   case 'T3'
    %%%%%%%%%%%%%%%% T3 THREE NODE TRIANGULAR ELEMENT %%%%%%%%%%%%%%%%%% 
    %   
    %               3
    %             /  \
    %            /    \
    %           /      \
    %          /        \
    %         /          \
    %        /            \
    %       /              \
    %      /                \
    %     /                  \
    %    1--------------------2
    %
    if size(coord,2) < 2
      disp('Error two coordinates needed for the T3 element')
    else
      xi=coord(1); eta=coord(2);
      N=[1-xi-eta;xi;eta];
      dNdxi=[-1,-1;1,0;0,1];
    end        

   case 'T6'
    %%%%%%%%%%%%%%%%%% T6 SIX NODE TRIANGULAR ELEMENT %%%%%%%%%%%%%%%%%%
    %   
    %               3
    %             /  \
    %            /    \
    %           /      \
    %          /        \
    %         6          5
    %        /            \
    %       /              \
    %      /                \
    %     /                  \
    %    1---------4----------2
    %
    if size(coord,2) < 2
      disp('Error two coordinates needed for the T6 element')
    else
      xi=coord(1); eta=coord(2);
      N=[1-3*(xi+eta)+4*xi*eta+2*(xi^2+eta^2);
                                  xi*(2*xi-1);
                                eta*(2*eta-1);
                              4*xi*(1-xi-eta);
                                     4*xi*eta;
                              4*eta*(1-xi-eta)];
        
      dNdxi=[4*(xi+eta)-3,   4*(xi+eta)-3;
                   4*xi-1,              0; 
                        0,        4*eta-1;
           4*(1-eta-2*xi),          -4*xi;
                    4*eta,           4*xi;
                   -4*eta,  4*(1-xi-2*eta)];
    end
    
    
   case 'Q4'
    %%%%%%%%%%%%%%% Q4 FOUR NODE QUADRILATERIAL ELEMENT %%%%%%%%%%%%%%%%
    %
    %    4--------------------3
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    1--------------------2
    %
    if size(coord,2) < 2
      disp('Error two coordinates needed for the Q4 element')
    else
      xi=coord(1); eta=coord(2);
      N=1/4*[ (1-xi)*(1-eta);
              (1+xi)*(1-eta);
              (1+xi)*(1+eta);
              (1-xi)*(1+eta)];
          
      dNdxi=1/4*[-(1-eta), -(1-xi);
		         1-eta,    -(1+xi);
		         1+eta,      1+xi;
                -(1+eta),   1-xi];
    end
    
    case 'Q8'
    %%%%%%%%%%%%%%% Q9 NINE NODE QUADRILATERIAL ELEMENT %%%%%%%%%%%%%%%%
    %
    %    4---------7----------3
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    8                    6
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    1----------5---------2
    %
    if size(coord,2) < 2
      disp('Error two coordinates needed for the Q8 element')
    else
      xi=coord(1); eta=coord(2);
      N=1/4*[-1*(1-xi)*(1-eta)*(1+xi+eta);
             -1*(1+xi)*(1-eta)*(1-xi+eta);
             -1*(1+xi)*(1+eta)*(1-xi-eta);
             -1*(1-xi)*(1+eta)*(1+xi-eta);
              2*(1-xi^2)*(1-eta);
              2*(1+xi)*(1-eta^2);
              2*(1-xi^2)*(1+eta);
              2*(1-xi)*(1-eta^2)];
             
      dNdxi=1/4*[(1-eta)*(2*xi+eta),(1-xi)*(2*eta+xi);
                 (1-eta)*(2*xi-eta),(1+xi)*(2*eta-xi);
                 (1+eta)*(2*xi+eta),(1+xi)*(2*eta+xi);
                 (1+eta)*(2*xi-eta),(1-xi)*(2*eta-xi);
                 -4*xi*(1-eta),-2*(1-xi^2);
                 2*(1-eta^2),-4*eta*(1+xi);
                -4*xi*(1+eta),2*(1-xi^2);
                -2*(1-eta^2),-4*eta*(1-xi)];                 
    end   
    
    
   case 'Q9'
    %%%%%%%%%%%%%%% Q9 NINE NODE QUADRILATERIAL ELEMENT %%%%%%%%%%%%%%%%
    %
    %    4---------7----------3
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    8          9         6
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    1----------5---------2
    %
    if size(coord,2) < 2
      disp('Error two coordinates needed for the Q9 element')
    else
      xi=coord(1); eta=coord(2);
      N=1/4*[xi*eta*(xi-1)*(eta-1);
             xi*eta*(xi+1)*(eta-1);
             xi*eta*(xi+1)*(eta+1);
             xi*eta*(xi-1)*(eta+1);
            -2*eta*(xi+1)*(xi-1)*(eta-1);
            -2*xi*(xi+1)*(eta+1)*(eta-1);
            -2*eta*(xi+1)*(xi-1)*(eta+1);
            -2*xi*(xi-1)*(eta+1)*(eta-1);
             4*(xi+1)*(xi-1)*(eta+1)*(eta-1)];
         
      dNdxi=1/4*[eta*(2*xi-1)*(eta-1),xi*(xi-1)*(2*eta-1);
                 eta*(2*xi+1)*(eta-1),xi*(xi+1)*(2*eta-1);
                 eta*(2*xi+1)*(eta+1),xi*(xi+1)*(2*eta+1);
                 eta*(2*xi-1)*(eta+1),xi*(xi-1)*(2*eta+1);
                 -4*xi*eta*(eta-1),-2*(xi+1)*(xi-1)*(2*eta-1);
                 -2*(2*xi+1)*(eta+1)*(eta-1),-4*xi*eta*(xi+1);
                 -4*xi*eta*(eta+1), -2*(xi+1)*(xi-1)*(2*eta+1);
                 -2*(2*xi-1)*(eta+1)*(eta-1),-4*xi*eta*(xi-1);
                  8*xi*(eta^2-1),      8*eta*(xi^2-1)];
              
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   otherwise
    disp(['Element ',type,' not yet supported'])
    N=[]; dNdxi=[];
  end
 
  I=eye(dim);
  Nv=[];
  for i=1:size(N,1)
    Nv=[Nv;I*N(i)];
  end
  
  if ( dim == 1 )
    B=dNdxi;
  elseif ( dim == 2 )
    B=zeros(dim*size(N,1),3);
    
    B(1:dim:dim*size(N,1)-1,1) = dNdxi(:,1);
    B(2:dim:dim*size(N,1),2)   = dNdxi(:,2);
    
    B(1:dim:dim*size(N,1)-1,3) = dNdxi(:,2);
    B(2:dim:dim*size(N,1),3)   = dNdxi(:,1);
  
    
  end
end    % end of function
  
 