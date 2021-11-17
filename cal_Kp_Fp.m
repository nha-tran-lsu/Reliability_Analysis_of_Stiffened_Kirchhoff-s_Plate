function [ Kp,Fp ] = cal_Kp_Fp(Load)
% global Load 
global sdof nel nnel 
global gcoord ele_nods H ndof edof

Kp=sparse(sdof,sdof);		% initialization of global stiffness matrix
Fp=sparse(sdof,1);

% 3-intergal point
point=[ 0.5 0.5 0;...
        0.5 0.0 0.5;...
        0.0 0.5 0.5];
weight=[1/3 1/3 1/3];
% Index for subelement
as=[2 3 4;...
    3 1 4;...
    1 2 4];
for iel=1:nel
   n1=ele_nods(iel,1);% extract nodes of element
   n2=ele_nods(iel,2);
   n3=ele_nods(iel,3);
   nd=[n1 n2 n3];
% Extract coorddinate of nodes
   x(1)=gcoord(n1,1);   y(1)=gcoord(n1,2); 
   x(2)=gcoord(n2,1);   y(2)=gcoord(n2,2);    
   x(3)=gcoord(n3,1);   y(3)=gcoord(n3,2); 
   
   loadnes=sparse(edof+3,1);        % subelement force vector
   stiff=sparse(edof+3,edof+3);     % subelement stiffness matrix 

   k1=sparse(edof,edof);     
   k2=sparse(edof,edof);
   k3=sparse(edof,edof);
% Routing for subelement
   for se=1:3  
     a=as(se,:);
     if (se==1)
        for igauss=1:3
          gauss=point(igauss,:);
          we=weight(igauss);  
          [ k1 ] = Stiff_HCT_Sub_T1_SP( edof,H,x,y,gauss,we,k1 )   ; 
        end
       f1=Load_HCT_Sub_T1_SP(Load,x,y);
       id=[a(1) a(2) a(3)];
       index=get_eledof(id,nnel,ndof);
       % assemble element-2 matrices 
       stiff(index,index)=stiff(index,index)+k1;
       loadnes(index)=loadnes(index)+f1;
     elseif (se==2)
        for igauss=1:3
          gauss=point(igauss,:);
          we=weight(igauss);
         [ k2 ] = Stiff_HCT_Sub_T2_SP( edof,H,x,y,gauss,we,k2 );    
        end
        f2=Load_HCT_Sub_T2_SP(Load,x,y);
        id=[a(1) a(2) a(3)];
        index=get_eledof(id,nnel,ndof);
        % assemble element-2 matrices 
        stiff(index,index)=stiff(index,index)+k2;
        loadnes(index)=loadnes(index)+f2;  
     else
        for igauss=1:3
          gauss=point(igauss,:);
          we=weight(igauss);
          [ k3 ] = Stiff_HCT_Sub_T3_SP( edof,H,x,y,gauss,we,k3 );        
        end
        f3=Load_HCT_Sub_T3_SP(Load,x,y);
        id=[a(1) a(2) a(3)];
        index=get_eledof(id,nnel,ndof);
        % assemble element-3 matrices 
        stiff(index,index)=stiff(index,index)+k3;
        loadnes(index)=loadnes(index)+f3;
     end
  end
 % eliminate node at centroid of triangle element
 C=eliminate_SP(x,y);
  kee=stiff(1:9,1:9);
  keo=stiff(1:9,10:12);
  koe=stiff(10:12,1:9);
  koo=stiff(10:12,10:12);
  fee=loadnes(1:9);
  foo=loadnes(10:12);
  ke=kee+keo*C+C'*koe+C'*koo*C;
  fe=fee+C'*foo;
  index=get_eledof(nd,nnel,ndof);
  % assemble element matrices 
  Kp(index,index)=Kp(index,index)+ke;
  Fp(index)=Fp(index)+fe;
end

end

