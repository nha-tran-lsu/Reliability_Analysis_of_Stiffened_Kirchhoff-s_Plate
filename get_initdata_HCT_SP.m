function [gcoord,ele_nods,bcdof,bcval]=get_initdata_HCT_SP
global Lx Hy nx ny
dx=Lx/nx;          % the length of side of elements in x-axis
dy=Hy/ny;          % the length of side of elements in y-axis
gcoord=[];
for j=1:ny+1
    for i=1:nx+1
        gcoord=[gcoord; (i-1)*dx (j-1)*dy;]; 
    end
end

for j=1:ny
    for i=1:nx
        e1=2*(j-1)*nx+i;
        e2=e1+nx;
% Node's position of e1
        ne1(1)=(j-1)*(nx+1)+i;
        ne1(2)=ne1(1)+1;
        ne1(3)=ne1(2)+nx+1;    
% Node's position of e2
        ne2(1)=ne1(1);
        ne2(2)=ne1(3);
        ne2(3)=ne2(2)-1;
% Matrix containing nodes of each element 
        ele_nods(e1,:)=[ne1(1) ne1(2) ne1(3)];
        ele_nods(e2,:)=[ne2(1) ne2(2) ne2(3)];
    end
end
% Boundary condition
L1 = find(gcoord(:,2)==min(gcoord(:,2)));  % at y = 0 (along X-axes);
L2 = find(gcoord(:,1)==max(gcoord(:,1)));  % at x = a (along Y-axes);
L3 = find(gcoord(:,2)==max(gcoord(:,2)));  % at y = b (along X-axes);
L4 = find(gcoord(:,1)==min(gcoord(:,1)));  % at x = 0 (along Y-axes);
n = length(L1) ;
m = length(L2) ;
% disp('plate is clamped at all the edges');
% dofL1 = zeros(1,3*n) ;
% dofL2 = zeros(1,3*m) ;
% dofL3 = zeros(1,3*n) ;
% dofL4 = zeros(1,3*m) ;
% for i = 1:n
%     i1 = 3*i ;      i2 = i1-1 ;     i3 = i1-2 ;
%     dofL1(i1) = 3*L1(i) ;       dofL1(i2) = 3*L1(i)-1 ;     dofL1(i3) = 3*L1(i)-2 ;
%     dofL3(i1) = 3*L3(i) ;       dofL3(i2) = 3*L3(i)-1 ;     dofL3(i3) = 3*L3(i)-2 ;
% end   
% for i = 1:m
%     i1 = 3*i ;      i2 = i1-1 ;     i3 = i1-2 ;
%     dofL2(i1) = 3*L2(i) ;       dofL2(i2) = 3*L2(i)-1 ;     dofL2(i3) = 3*L2(i)-2 ;
%     dofL4(i1) = 3*L4(i) ;       dofL4(i2) = 3*L4(i)-1 ;     dofL4(i3) = 3*L4(i)-2 ;
% end
%     L1UL3 = union(dofL1,dofL3) ; 
%     L2UL4 = union(dofL2,dofL4) ;
%     bcdof = union(L1UL3,L2UL4) ;
% %     bcval = 0;
%     bcval = zeros(length(bcdof),1) ;
% end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % 
% disp('plate is simply Supported at all the edges');
% Qui uoc , w_x=-dw/dy, w_y=dw/dx
dofL1 = zeros(1,2*n) ;
dofL2 = zeros(1,2*m) ;
dofL3 = zeros(1,2*n) ;
dofL4 = zeros(1,2*m) ;
for i = 1:n
    i1 = 2*i ;      i2 = i1-1 ;
    dofL1(i1) = 3*L1(i) ;       dofL1(i2) = 3*L1(i)-2 ;
    dofL3(i1) = 3*L3(i) ;       dofL3(i2) = 3*L3(i)-2 ;
end   
for i = 1:m
    i1 = 2*i ;      i2 = i1-1 ;
    dofL2(i1) = 3*L2(i)-1 ;       dofL2(i2) = 3*L2(i)-2 ;
    dofL4(i1) = 3*L4(i)-1 ;       dofL4(i2) = 3*L4(i)-2 ;
end
    L1UL3 = union(dofL1,dofL3) ; 
    L2UL4 = union(dofL2,dofL4) ;
    bcdof = union(L1UL3,L2UL4) ;
%     bcval = 0;
    bcval = zeros(length(bcdof),1) ;
end