function [CP] = plastic_mat(E,nu,phi,tsi,iel,kk,stress)

% Generates the plastic constitutive matix CP at each gauss point

        s_xx=stress(1,kk,iel);s_yy=stress(2,kk,iel);
        t_xy=stress(3,kk,iel);s_zz=stress(4,kk,iel);
        p1=(s_xx+s_yy+s_zz)/3;

        t=sqrt((s_xx-s_yy)^2+(s_yy-s_zz)^2+(s_zz-s_xx)^2+6*t_xy^2)/sqrt(3);
        q=sqrt(1.5)*t;
        sx=s_xx-p1; 
        sy=s_yy-p1;
        sz=s_zz-p1;
        j2=-sx*sy-sy*sz-sz*sx+t_xy^2;
        j3=sx*sy*sz-(sx*t_xy^2);
        th=-3*sqrt(3)*j3/(2*j2^1.5);
                 if th>1
                    th=1;
                 end
                 if th<-1
                    th=-1;
                 end
           theta1=asind(th)/3;
  
if sind(theta1)>0.495         % close to theta=30 corner,smoothen the curve 
                              % with triaxial compression case s1 > s2 = s3
    sw=-1; 
    cphi=0.25*sqrt(3/j2)*(1+(sw*sind(phi)/3));
    ctsi=0.25*sqrt(3/j2)*(1+(sw*sind(tsi)/3));
    kphi=sind(phi)*(1+nu)/3;
    ktsi=sind(tsi)*(1+nu)/3;
    c1=kphi+cphi*(sx*(1-nu)+nu*(sy+sz));
    c2=kphi+cphi*(sy*(1-nu)+nu*(sx+sz));
    c3=cphi*(1-2*nu)*t_xy;
    c4=kphi+cphi*(sz*(1-nu)+nu*(sx+sy));
    r1=ktsi+ctsi*(sx*(1-nu)+nu*(sy+sz));
    r2=ktsi+ctsi*(sy*(1-nu)+nu*(sx+sz));
    r3=ctsi*(1-2*nu)*t_xy;
    r4=ktsi+ctsi*(sz*(1-nu)+nu*(sx+sy));
    A=[r1*c1 r1*c2 r1*c3 r1*c4;...
       r2*c1 r2*c2 r2*c3 r2*c4;... 
       r3*c1 r3*c2 r3*c3 r3*c4;...
       r4*c1 r4*c2 r4*c3 r4*c4];
   CP=E*A/((1+nu)*(1-2*nu)*(kphi*sind(tsi)+2*cphi*ctsi*j2*(1-2*nu)));
  
elseif -1*sind(theta1)> 0.49
     sw=1;                   % close to theta=-30 corner,smoothen the curve
                            % with triaxial extension case s1 = s2 > s3
    cphi=0.25*sqrt(3/j2)*(1+(sw*sind(phi)/3));
    ctsi=0.25*sqrt(3/j2)*(1+(sw*sind(tsi)/3));
    kphi=sind(phi)*(1+nu)/3;
    ktsi=sind(tsi)*(1+nu)/3;
    c1=kphi+cphi*(sx*(1-nu)+nu*(sy+sz));
    c2=kphi+cphi*(sy*(1-nu)+nu*(sx+sz));
    c3=cphi*(1-2*nu)*t_xy;
    c4=kphi+cphi*(sz*(1-nu)+nu*(sx+sy));
    r1=ktsi+ctsi*(sx*(1-nu)+nu*(sy+sz));
    r2=ktsi+ctsi*(sy*(1-nu)+nu*(sx+sz));
    r3=ctsi*(1-2*nu)*t_xy;
    r4=ktsi+ctsi*(sz*(1-nu)+nu*(sx+sy));
    A=[r1*c1 r1*c2 r1*c3 r1*c4;...
       r2*c1 r2*c2 r2*c3 r2*c4;... 
       r3*c1 r3*c2 r3*c3 r3*c4;...
       r4*c1 r4*c2 r4*c3 r4*c4];
   CP=E*A/((1+nu)*(1-2*nu)*(kphi*sind(tsi)+2*cphi*ctsi*j2*(1-2*nu)));
   
else                                               % all other cases
    alpha=atand(abs((s_xx-s_yy)/(2*t_xy)));
    k1=1;k2=1;
        if abs(s_xx)>abs(s_yy)
            k1=-1;
        end
        if t_xy <0
            k2=-1;
        end
    c1=sind(phi)+k1*(1-2*nu)*sind(alpha);
    c2=sind(phi)-k1*(1-2*nu)*sind(alpha);
    c3=k2*(1-2*nu)*cosd(alpha);
    c4=2*nu*sind(phi);
    r1=sind(tsi)+k1*(1-2*nu)*sind(alpha);
    r2=sind(tsi)-k1*(1-2*nu)*sind(alpha);
    r3=k2*(1-2*nu)*cosd(alpha);
    r4=2*nu*sind(tsi);
    A=[r1*c1 r1*c2 r1*c3 r1*c4;...
       r2*c1 r2*c2 r2*c3 r2*c4;... 
       r3*c1 r3*c2 r3*c3 r3*c4;...
       r4*c1 r4*c2 r4*c3 r4*c4];
   CP=E*A/(2*(1+nu)*(1-2*nu)*(1-2*nu+sind(phi)*sind(tsi)));
   
end
end    % end of function