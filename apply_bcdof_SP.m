function [GK GP]=apply_bcdof_SP(GK,GP,bcdof,bcval)
n=length(bcdof);
 sdof=size(GK);
 for i=1:n
    c=bcdof(i);
    for j=1:sdof
    GK(c,j)=0;
    end
    GK(c,c)=1;
%     GP(c) = bcval;
    GP(c)=bcval(i);
    
 end
end

