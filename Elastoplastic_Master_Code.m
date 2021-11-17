%  This script file is a FE code for a Rectangual element with 
%  4,8 or 9 nodes or for a triangular element with 3 or 6 nodes 
%  to analyze stability of slopes  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
tic; 
feval('setpath')
global node element elemType 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E=60e6;
nu=0.25;
L= 50;
D= 30;

phi=30;
tsi=30;
c=1e3;
numx=10;
numy=20;
s_v=0;
s_h=0;  
elemType = 'Q4' ;
normal_order =2; 
df=0;
nsteps=100; 
maxit=20;
tol=1e-2;  
load_edge1=0;  % set to be zero
load_edge2=50; %L;  % set to be the whole length
%%%% New line 
  gamma = [10 20  20;
           20 30 25;
          30 40 28];
    % each row has values:: y1 y2  gamma
    % where y = y1 represents the lower boundary line of segment 
    %       y = y2 represents the upper boundary line of segment
%%%%%
%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C= E/(1+nu)/(1-2*nu)*[ 1-nu nu 0 nu;nu 1-nu 0 nu;0 0 0.5-nu 0;nu nu 0 1-nu]; 
dt=4*(1+nu)*(1-2*nu)/(E*(1-2*nu+(sind(phi))^2));  % pseudo time stepping parameter
% Four corner points
pt1 = [0 10];
pt2 = [50  10];
pt3 = [ 50   18];
%%%% Additional
%%%% Conrner points at top layer matrix pte1 is introduced
%pte1= [50 210;  71  180; 75 180; 95 155;105  155; 155  105];% corner points at topside
pte1 = [10 40; 25 25; 35  25; 45 20; 50 18 ];
%%
pt4 = [ 0  40];

%%%% Additional :: Plotting of surface after addition of corner pts
xkk=[pt1(:,1)  pt4(:,1)  [pte1(:,1)]'  pt3(:,1)   pt2(:,1)  pt1(:,1) ];
ykk=[pt1(:,2)  pt4(:,2)  [pte1(:,2)]'  pt3(:,2)   pt2(:,2)  pt1(:,2)];
plot([xkk],[ykk],'LineWidth',3);
hold on;
%%%%%%%

sigmatox=0;            % make the horizontal load zero
disp([num2str(toc),'    MESH GENERATION....'])
 switch elemType
     case {'Q9','Q8','Q4','T3'}
      [node,element] = mesh_region(pt1, pt2, pt3, pte1, pt4,numx,numy,elemType);
     otherwise 
      [node,element] =mesh_t6_elem(L,D,numx,numy);
 end
[topEdge,topEdge1,dispNodes,dispNodes1,leftNodes1]=supportcond(elemType,numx,numy); 
dispNodes2=dispNodes(2:end);
dispNodes=dispNodes(1);
dispNodes1=leftNodes1;
numnode = size(node,1);
numelem = size(element,1);
nonelm=size(element,2);
total_unknown = 2*numnode;  
udofs  = [dispNodes;(dispNodes1.*2)-1]; 
vdofs  = [dispNodes.*2;dispNodes2.*2];                         
dofs=union(udofs(:),vdofs(:));                 
unknowndof=setdiff((1:total_unknown)',dofs);

disp([num2str(toc),'    STIFFNESS MATRIX COMPUTATION....'])
K=stiffness_matrix(node,element,elemType,normal_order,C);

stress=zeros(4,nonelm,numelem); % generate initial stress
stress(1,:,:)=s_v;
stress(2,:,:)=s_h;
stress(3,:,:)=0;
stress(4,:,:)=nu*(s_v+s_h);
[p,q,theta] =invariants1(element,stress);
load=s_v;

strain=zeros(4,nonelm,numelem); % set variables to zero
strainP=zeros(4,nonelm,numelem);
stress_tr=zeros(4,nonelm,numelem);
ui=zeros(total_unknown,1);
force=zeros(total_unknown,1);

f_old=zeros(total_unknown,1);

r=zeros(total_unknown,1);
b=zeros(total_unknown,1);
du=zeros(total_unknown,1);
dgds=zeros(4,nonelm,numelem);

%%%%% Introduce the new line  ::: Selfweight is initial force
%%% The value of force get updated by selfweight

selfwt=selfwt_matrix(elemType,normal_order,gamma,node,element);

%%%%%%%%%%%%%%%%%
                              % prepare space for plotting data
pvq=zeros(2,nsteps+1);     pvq(1,1)=p(1);    pvq(2,1)=q(1);
evsyy=zeros(2,nsteps+1);   evsyy(1,1)=0;     evsyy(2,1)=stress(2,1,1);
epsvq=zeros(2,nsteps+1);   epsvq(1,1)=0;     epsvq(2,1)=q(1);
fvu=zeros(2,nsteps+1);     fvu(1,1)=0;       fvu(2,1)=load;   

 kk1 = 1;
 
 % start load stepping
for steps=1:nsteps
    stepno=steps;
    err=1; nit=0;  
    sigmatoy=df;
    switch elemType
     case {'Q9','Q8','T6'}          
          [f,sctry]=force_matrix689(node,topEdge,topEdge1,sigmatoy,sigmatox,load_edge1,load_edge2);
     case{'Q4','T3'}
          [f,sctry]=force_matrix(node,topEdge,sigmatoy,sigmatox,load_edge1,load_edge2);  
    end  
    
%%%%%%
F1=f(unknowndof)+selfwt(unknowndof);
%%%%%%
    r=zeros(total_unknown,1);    % reset unknown variables to zero
    DEPS_PLA=zeros(4,nonelm,numelem);
    Du=zeros(total_unknown,1);
    du_old=zeros(total_unknown,1);
    Deps=zeros(4,nonelm,numelem);
    dsig_pla=zeros(4,nonelm,numelem);
    Deps_ela=zeros(4,nonelm,numelem);
    Dsig=zeros(4,nonelm,numelem);
    deps_pla=zeros(4,nonelm,numelem);
   
                     % start iteration loop
while (err>tol) && (nit<maxit)
        nit=nit+1;
% Modification done::du(unknowndof)=K(unknowndof,unknowndof)\f(unknowndof); 
         du(unknowndof)=K(unknowndof,unknowndof)\F1; 
        Du=Du+du;
for iel = 1 : numelem       % start looping on all elements
    sctr = element(iel,:); % element connectivity
    nn   = length(sctr);   % number of nodes per element
    eldof =elementdof(elemType,sctr);  
   [W,Q] = gauss_rule(iel,elemType,normal_order);     
       
    for kk = 1:nn       % start looping on all gauss points
        pt = Q(kk,:);                             
        [N,dNdxi] = shape_func(elemType,pt);  
        J0 = node(sctr,:)'*dNdxi;                    
        Bfem4 =Bmatrix4(pt,elemType,iel);
        Deps(:,kk,iel)=Bfem4*du(eldof);      
        Deps(end,:)=0;
        Deps_ela(:,kk,iel)=Deps(:,kk,iel)-DEPS_PLA(:,kk,iel);
        Dsig(:,kk,iel)= C*Deps_ela(:,kk,iel);
        stress_tr(:,kk,iel)=stress(:,kk,iel)+Dsig(:,kk,iel);        
        [p,q,theta] =invariants2(kk,iel,stress_tr);
        F=p*sind(phi)+q*((cosd(theta)/sqrt(3))-(sind(theta)*sind(phi)/3))-c*cosd(phi);
      
            if F<0  
               stress(:,kk,iel)=stress_tr(:,kk,iel);
               err=0;
            else      
               [m1,m2,m3] = formm(kk,iel,stress_tr);
               [dg1,dg2,dg3] =formdg(tsi,q,theta);
               dgds(:,kk,iel)=(dg1*m1+dg2*m2+dg3*m3)*stress_tr(:,kk,iel);
               deps_pla(:,kk,iel)=dt*F*dgds(:,kk,iel);        
               DEPS_PLA(:,kk,iel)=DEPS_PLA(:,kk,iel)+deps_pla(:,kk,iel); 
               dsig_pla(:,kk,iel)=C*deps_pla(:,kk,iel);   
             
                       if nit==1
                          err=1;
                       else
                          err=max(abs(du_old(2:2:end)-du(2:2:end)));                           
                       end                     
            end  
           r(eldof) =r(eldof)+Bfem4'*(C*deps_pla(:,kk,iel))*W(kk)*det(J0);  
    end                 % end of looping on GPs           
          f(sctry)= f(sctry)+r(sctry);  
end                     % end of looping on elements 
         
                           
         du_old=du;
end                   %end of iteration

     stress(:,kk,iel)=stress(:,kk,iel);
     s=permute(stress,[2 1 3]);  sigma=permute(stress,[3 2 1]);   
     ui=ui+Du;
     u_x=ui(1:2:2*numnode-1) ;
     u_y=ui(2:2:2*numnode) ;
     strain=strain+Deps;   
     strainP=strainP+DEPS_PLA;        
     ule=numelem-numx+1;
     xxx=strain(2,1,ule);
     yyy=stress(2,1,ule);               
     eps_vol=strainP(1,1,ule)+strainP(2,1,ule)+strainP(4,1,ule);
     stressule=stress(:,1,ule);   
     [p,q,theta] = invariants(stressule);
     pvq(1,steps+1)=p(end,:);
     pvq(2,steps+1)=q(end,:);
     epsvq(1,steps+1)=eps_vol(end,:);
     epsvq(2,steps+1)=q(end,:);    
     evsyy(1,steps+1)=xxx(end,:);
     evsyy(2,steps+1)=yyy(end,:);
     load=load+df;        
     fvu(1,steps+1)=u_y(end,:);
     fvu(2,steps+1)=load(end,:);
end                     %end of load step

figure
hold on
plot(abs(pvq(1,:)),pvq(2,:),'--b*','linewidth',2);

  xlabel({'P'},'FontSize',16);
  ylabel({'q'},'FontSize',16);

figure
hold on
plot(abs(epsvq(1,:)),epsvq(2,:),'--r*','linewidth',2);

  xlabel({'epsvol'},'FontSize',16);
  ylabel({'q'},'FontSize',16);

  
  figure
hold on
plot(abs(evsyy(1,:)),abs(evsyy(2,:)),'--r*','linewidth',2);

  xlabel({'epsyy'},'FontSize',16);
  ylabel({'sigyy'},'FontSize',16);
  
  figure
hold on
plot(abs(fvu(1,:)),abs(fvu(2,:)),'--r*','linewidth',2);

  xlabel({'u_y'},'FontSize',16);
  ylabel({'load'},'FontSize',16);
  
% Plot the FEM mesh 
dispNodes1=union(dispNodes1,dispNodes2);
plot_m(elemType,dispNodes,dispNodes1)
title('Undeformed FE mesh')
%Plot numerical deformed configuration
dispnorm=L/max(sqrt(u_x.^2+u_y.^2));
fac =dispnorm*0.05;                    %magnification factor
plot_def(fac,u_x,u_y,elemType,dispNodes,dispNodes1);
%plot stress and deformation intensity with a colormap
plot_defo(fac,u_x,u_y,elemType)
plot_sig(fac,u_x,u_y,elemType,sigma)