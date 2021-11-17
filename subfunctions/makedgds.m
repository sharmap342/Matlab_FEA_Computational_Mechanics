function [dgds] =makedgds(tsi,element,iel,kk,stress)

%Forms the partial derivatives of the plastic potential function with
%respect to stress at each gauss point.

numelem=size(element,1);
nonelm=size(element,2);
dgds=zeros(4,nonelm,numelem);

        s_xx=stress(1,kk,iel);s_yy=stress(2,kk,iel);
        t_xy=stress(3,kk,iel);s_zz=stress(4,kk,iel);
        p1=(s_xx+s_yy+s_zz)/3;

        t=sqrt((s_xx-s_yy)^2+(s_yy-s_zz)^2+(s_zz-s_xx)^2+6*t_xy^2)/sqrt(3);
        q1=sqrt(1.5)*t;
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
  
 if sind(theta1)>0.49  % close to theta=30 corner,smoothen the curve
                       % with triaxial compression case s1 > s2 = s3
    sw=-1;     
    ctsi=0.25*sqrt(3/j2)*(1+(sw*sind(tsi)/3));  
    
    dgds(:,kk,iel)=[(sind(tsi)/3)+sx*ctsi;...
                    (sind(tsi)/3)+sy*ctsi;...
                    t_xy*ctsi;...
                    (sind(tsi)/3)+sz*ctsi];
        
 elseif (-1*sind(theta1))>0.49  % close to theta=-30 corner,smoothen the curve
                                % curve with triaxial extension case s1 = s2 > s3
    sw=1;   
    ctsi=0.25*sqrt(3/j2)*(1+(sw*sind(tsi)/3));
         
    dgds(:,kk,iel)=[(sind(tsi)/3)+sx*ctsi;...
                    (sind(tsi)/3)+sy*ctsi;...
                    t_xy*ctsi;...
                    (sind(tsi)/3)+sz*ctsi];
 else                            
    alpha=atand(abs((s_xx-s_yy)/(2*t_xy)));  % all other conditions
    k1=1;k2=1;
                     if abs(s_xx)>abs(s_yy)
                         k1=-1;
                     end
                     if t_xy <0
                         k2=-1;
                     end
   
    dgds(:,kk,iel)=[sind(tsi)+k1*sind(alpha);...
                    sind(tsi)-k1*sind(alpha);...
                    2*k2*cosd(alpha);...
                             0];
 end
 
end % end of function