function [ Ksy ] = cal_Ksy(Emodule)
global Hy ny ISx
% global Emodule
nnelSy=2;               % Number of nodes per element
ndofSy=2;               % Number of dofs per node
nnodeSy=(ny+1);         % Total number of nodes in system
sdofSy=nnodeSy*2;         % Total degrees of freedom in system
dy=Hy/ny;          % the length of side of elements in y-axis
ycoordS=0:dy:Hy;
for i=1:ny
    nodesSy(i,:)=[i i+1];
end
Ksy=sparse(sdofSy,sdofSy);	
for iel=1:ny
   n1=nodesSy(iel,1);      n2=nodesSy(iel,2); % extract nodes of element
   nd=[n1 n2];
   % Extract coorddinate of nodes
   y(1)=ycoordS(n1);    y(2)=ycoordS(n2);
   [ Kesy ] = Stiff_Beam_TVN108_oy( y,Emodule,ISx);
   index=get_eledof(nd,nnelSy,ndofSy);
   Ksy(index,index)=Ksy(index,index) + Kesy;
end
end

