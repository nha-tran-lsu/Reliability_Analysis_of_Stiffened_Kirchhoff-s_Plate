function [ Tsx, Tsy ] = Assembling_Ksx_Ksy
global gcoord sdof
Sx = find(gcoord(:,2)==max(gcoord(:,2))/2);  % at y = Hy/2 (along X-axes);
Sy = find(gcoord(:,1)==max(gcoord(:,1))/2);  % at x = Lx/2 (along Y-axes);
n = length(Sx) ;
m = length(Sy) ;
dofSx = zeros(1,2*n) ;
dofSy = zeros(1,2*m) ;
for i = 1:n
    i1 = 2*i ;      i2 = i1-1 ;
    dofSx(i1)=3*Sx(i);    dofSx(i2)=3*Sx(i)-2; 
%     dw/dx                     w
end   
Tsx = sparse(2*n,sdof);
for i=1:2*n
    Tsx(i,dofSx(i))=1;
end
for j=1:m
    j1 = 2*j ;      j2 = j1-1 ;
    dofSy(j1)=3*Sy(j)-1;   dofSy(j2)=3*Sy(j)-2;
%     -dw/dy                     w
end
Tsy = sparse(2*m,sdof);
for j=1:2*m
    Tsy(j,dofSy(j))=1;
end
end

