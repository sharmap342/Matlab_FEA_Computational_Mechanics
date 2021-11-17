function [m1,m2,m3] = formm(kk,iel,stress)

%Generates the partial derivatives of the p, J2 and J3 with respect to
%stress for use in the plastic potential derivative.

  s_xx=stress(1,kk,iel);s_yy=stress(2,kk,iel);
  t_xy=stress(3,kk,iel);s_zz=stress(4,kk,iel);
  p1=(s_xx+s_yy+s_zz)/3;

  sx=s_xx-p1; 
  sy=s_yy-p1;
  sz=s_zz-p1;
  
  m1=1/(3*(s_xx+s_yy+s_zz))*[1 1 0 1;1 1 0 1;0 0 0 0;1 1 0 1];
  
  m2=1/3*[2 -1 0 -1;-1 2 0 -1;0 0 6 0;-1 -1 0 2];
  
  m3=1/3*[sx sz t_xy sz;sz sy t_xy sx;t_xy t_xy -3*sz -2*t_xy;sz sx -2*t_xy sz];

end % emd of function