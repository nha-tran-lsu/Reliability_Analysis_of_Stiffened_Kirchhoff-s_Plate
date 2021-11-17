function [ Ksx ] = cal_Ksx(Emodule)
global Lx nx ISy
% global Emodule
nnelSx=2;               % Number of nodes per element
ndofSx=2;               % Number of dofs per node
nnodeSx=(nx+1);         % Total number of nodes in system
sdofSx=nnodeSx*2;         % Total degrees of freedom in system
dx=Lx/nx;          % the length of side of elements in x-axis
xcoordS=0:dx:Lx;
for i=1:nx
    nodesSx(i,:)=[i i+1];
end
Ksx=sparse(sdofSx,sdofSx);	
for iel=1:nx
   n1=nodesSx(iel,1);      n2=nodesSx(iel,2); % extract nodes of element
   nd=[n1 n2];
   % Extract coorddinate of nodes
   x(1)=xcoordS(n1);    x(2)=xcoordS(n2);
   [ Kesx ] = Stiff_Beam_TVN108_ox( x,Emodule,ISy);
   index=get_eledof(nd,nnelSx,ndofSx);
   Ksx(index,index)=Ksx(index,index) + Kesx;
end
end

