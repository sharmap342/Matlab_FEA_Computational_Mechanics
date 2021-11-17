function [p,q,theta] =invariants1(element,stress)

% calculates the stress invariants at each gauss point.This can be used  
% outside the iteration loop.

numelem=size(element,1);
nonelm=size(element,2);
p=zeros(1,nonelm,numelem);
q=zeros(1,nonelm,numelem);
theta=zeros(1,nonelm,numelem);

for iel = 1 : numelem 
    sctr = element(iel,:); % element connectivity
    nn   = length(sctr);   % number of nodes per element   
       
    for kk = 1:nn       
        s_xx=stress(1,kk,iel);s_yy=stress(2,kk,iel);
        t_xy=stress(3,kk,iel);s_zz=stress(4,kk,iel);
        p1=(s_xx+s_yy+s_zz)/3;
        t=sqrt((s_xx-s_yy)^2+(s_yy-s_zz)^2+(s_zz-s_xx)^2+6*t_xy^2)/sqrt(3);
        q1=sqrt(1.5)*t;
        sx=s_xx-p1; 
        sy=s_yy-p1;
        sz=s_zz-p1;
               if q1<1e-6
                  theta11=0;
               else   
                  j3=sx*sy*sz-(sx*t_xy^2);
                  sine=-3*j3*sqrt(6)/t^3;
                          if sine>1;
                             sine=1;
                          end
                          if sine<-1;
                             sine=-1;
                          end
                   theta11=1/3*(asind(sine));
               end
      p(:,kk,iel)=p1;
      q(:,kk,iel)=q1;
      theta(:,kk,iel)=theta11;              
    end                 % end of looping on GPs                             
end                     % end of looping on elements
end                     % end of function