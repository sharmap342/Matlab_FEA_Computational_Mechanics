function node=square_node_array(pt1,pt2,pt3,pte1,pt4,nnx,nny)

% Generates a quadratleral array of nodes between the counterclockwise 
% ordering of nodes pt1 - pt4,given number of elements in x and y direction 


if ( nargin < 6 )
   disp('Not enough parameters specified for quare_node_array function')

end

% get node spacing along u direction

  xi_pts=linspace(-1,1,nnx);
%   dx= xi_pts(2)-xi_pts(1);
%   xi_pts(nnx-1) = xi_pts(nnx-1)+ 0.95*dx;
%   xi_pts(nnx-2) = xi_pts(nnx-2)+ 0.8*dx;
%   xi_pts(nnx-3) = xi_pts(nnx-3)+ 0.6*dx;
%   xi_pts(nnx-3) = xi_pts(nnx-3)+ 0.8*dx;
% get node spacing along v direction

  eta_pts=linspace(-1,1,nny);


pt1m= pt1; % Do not change Lenght must be fixed
pt2m= pt2; % Do not change Lenght must be fixed
pt3m=pt3; % Top Layer corner need to modified as per slope
pt4m=pt4; % Top Layer corner need to modified as per slope

x_pts=[pt1m(1),pt2m(1),pt3m(1),pt4m(1)];
y_pts=[pt1m(2),pt2m(2),pt3m(2),pt4m(2)];

for k = 1: size(pte1,1) 
xe1(1,k) = -1 + 2*(pte1(k,1)-pt4(1))/(pt3(1)-pt4(1));
% isoparametric coordinate of pte1
end
    
TopX = [pt4(1); pte1(:,1); pt3(1)];
TopY =  [pt4(2); pte1(:,2); pt3(2)];
Tempxei = [-1 xe1 1];

for r=1:nny
  eta=eta_pts(r);
  k1 =1;
  for c = 1:nnx
    xi=xi_pts(c);
        % get interpolation basis at xi, eta
    N=shape_func('Q4',[xi,eta]);
    N=N(:,1);
        %%%%% additional code
        pt1m = pt1;
        pt2m = pt2; 
        % if xi>=0 && xi<=1
        pt3m(1) = pt3(1);
        pt4m(1) = pt4(1);
  
        if xi >= Tempxei(k1) && xi < Tempxei(k1+1)   
        %  Need to modify pt3 and pt4 as per corner points on Top Side
        else
        k1 = k1+1; 
        end
        % Need to modify pt3 and pt4 as per corner points on Top Side
         if xi == -1
                pt3m(2) = pt3(2);
                pt4m(2) = pt4(2);
        elseif xi == 1
               pt3m(2) = pt3(2);
               pt4m(2) = pt4(2);
         else
                m1 =( TopY(k1+1)- TopY(k1)) / ( TopX(k1+1)- TopX(k1));
                c1 =( TopX(k1+1)*TopY(k1) - TopY(k1+1)*TopX(k1) ) / ...
                   ( TopX(k1+1)- TopX(k1));  
                pt3m(2) = m1*pt3m(1) + c1;
                pt4m(2) = m1*pt4m(1) + c1;
        
        end
    
    
    x_pts=[pt1m(1),pt2m(1),pt3m(1),pt4m(1)];
    y_pts=[pt1m(2),pt2m(2),pt3m(2),pt4m(2)];
    %%%%
    node((r-1)*nnx+c,:)=[x_pts*N,y_pts*N];
  end
  
end
end

