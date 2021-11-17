function[tam]=Showmesh_Triangle(coord,ele)
fh = figure ;
set(fh,'name','Preprocessing for FEA','numbertitle','off');%,'color','w') ;
N=size(coord,1); Ne=max(size(ele)); tam=zeros(Ne,2); %hold on
for i=1:Ne
    if iscell(ele)==0
       node=ele(i,:); node(find(node==0))=[]; nnel=length(node);
    else 
       node=ele{i}; nnel=length(ele{i});
    end 
    x=coord(node,1);    y=coord(node,2);
    xx=sum(x,1)/nnel;      yy=sum(y,1)/nnel;    tam(i,:)=[xx yy];
    patch(x,y,'w');     axis('equal');
    text(xx,yy,num2str(i),'Color',[.8 0 0]);
end
text(coord(1:N,1),coord(1:N,2),num2str([1:N]'),'Color',[0 0 0.8]);

for i=1:length(tam)
    hold on
    plot(tam(i,1),tam(i,2),'linewidth',2);
end
title('Finite Element Mesh of Plate') ;
hold on
axis image