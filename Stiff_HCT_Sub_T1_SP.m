function [ ke ] = Stiff_HCT_Sub_T1_SP( edof,H,x,y,gauss,we,ke )
L1=gauss(1);    L2=gauss(2);    L3=gauss(3);
x1=x(1);    x2=x(2);    x3=x(3);    y1=y(1);    y2=y(2);    y3=y(3);
x0=(x1+x2+x3)/3;    y0=(y1+y2+y3)/3;
% A12=2*area_Sub_triangle_1
% o: The centroid of triangle element

%                                          3-node (2)
%                                            o
%                                           / \

%								L1=0      /		\   L3=0

%										/ Sub-T_1 \

%                         0-node  (3) o-------------o 2-node (1)
%                                            L2=0 
b30=y3-y0;                  b02=y0-y2;                  b23=y2-y3;
c30=x3-x0;                  c02=x0-x2;                  c23=x2-x3;
l30=sqrt(b30^2+c30^2);      l02=sqrt(b02^2+c02^2);      l23=sqrt(b23^2+c23^2);
% mu30=(l23^2-l02^2)/l30^2;   mu02=(l30^2-l23^2)/l02^2;   
mu23=(l02^2-l30^2)/l23^2;
A12=b30*(-c02)-(-c30)*b02;      % A12=2*area_Sub_triangle_1
A14=2*A12;


Q=[   b30   b02    b23 ;...
     -c30  -c02   -c23];
% Qui uoc , w_x=-dw/dy, w_y=dw/dx
 
D2L11 =[ 6*L1 - (6*L3*(b02*c30 - b30*c02))/(b02*c23 - b23*c02) + (6*L2*(b23*c30 - b30*c23))/(b02*c23 - b23*c02), (2*A12*L3*b02)/(b02*c23 - b23*c02) - (2*A12*L2*b23)/(b02*c23 - b23*c02), (2*A12*L2*c23)/(b02*c23 - b23*c02) - (2*A12*L3*c02)/(b02*c23 - b23*c02), 0, 0, 0, 0, 0, 0];
D2L22 =[ 0, 0, 0, 6*L2 + (6*L1*(b02*c23 - b23*c02))/(b23*c30 - b30*c23) - (6*L3*(b02*c30 - b30*c02))/(b23*c30 - b30*c23), (2*A12*L1*b23)/(b23*c30 - b30*c23) - (2*A12*L3*b30)/(b23*c30 - b30*c23), (2*A12*L3*c30)/(b23*c30 - b30*c23) - (2*A12*L1*c23)/(b23*c30 - b30*c23), 0, 0, 0];
D2L33 =[ 0, 0, 0, 0, 0, 0, 6*L3 - (6*L1*(b02*c23 - b23*c02))/(b02*c30 - b30*c02) - (6*L2*(b23*c30 - b30*c23))/(b02*c30 - b30*c02), (2*A12*L1*b02)/(b02*c30 - b30*c02) - (2*A12*L2*b30)/(b02*c30 - b30*c02), (2*A12*L2*c30)/(b02*c30 - b30*c02) - (2*A12*L1*c02)/(b02*c30 - b30*c02)];
D2L12 =[ (3*L3*(b02*c23 - b23*c02 + 2*b02*c30 - 2*b30*c02 + 3*b23*c30 - 3*b30*c23 - b02*c23*mu23 + b23*c02*mu23 - b23*c30*mu23 + b30*c23*mu23))/(2*(b02*c23 - b23*c02)) + (6*L1*(b23*c30 - b30*c23))/(b02*c23 - b23*c02), (A14*L3*c23)/l23^2 - (A12*L3*(2*b02 + 3*b23 - b23*mu23))/(2*(b02*c23 - b23*c02)) - (2*A12*L1*b23)/(b02*c23 - b23*c02), (A12*L3*(2*c02 + 3*c23 - c23*mu23))/(2*(b02*c23 - b23*c02)) + (A14*L3*b23)/l23^2 + (2*A12*L1*c23)/(b02*c23 - b23*c02), (3*L3*(3*b02*c23 - 3*b23*c02 + 2*b02*c30 - 2*b30*c02 + b23*c30 - b30*c23 + b02*c23*mu23 - b23*c02*mu23 + b23*c30*mu23 - b30*c23*mu23))/(2*(b23*c30 - b30*c23)) + (6*L2*(b02*c23 - b23*c02))/(b23*c30 - b30*c23), (A12*L3*(3*b23 + 2*b30 + b23*mu23))/(2*(b23*c30 - b30*c23)) + (A14*L3*c23)/l23^2 + (2*A12*L2*b23)/(b23*c30 - b30*c23), (A14*L3*b23)/l23^2 - (A12*L3*(3*c23 + 2*c30 + c23*mu23))/(2*(b23*c30 - b30*c23)) - (2*A12*L2*c23)/(b23*c30 - b30*c23), 0, 0, 0];
D2L13 =[ (3*L2*(b02*c23 - b23*c02 + 2*b02*c30 - 2*b30*c02 + 3*b23*c30 - 3*b30*c23 - b02*c23*mu23 + b23*c02*mu23 - b23*c30*mu23 + b30*c23*mu23))/(2*(b02*c23 - b23*c02)) - (6*L1*(b02*c30 - b30*c02))/(b02*c23 - b23*c02), (A14*L2*c23)/l23^2 - (A12*L2*(2*b02 + 3*b23 - b23*mu23))/(2*(b02*c23 - b23*c02)) + (2*A12*L1*b02)/(b02*c23 - b23*c02), (A12*L2*(2*c02 + 3*c23 - c23*mu23))/(2*(b02*c23 - b23*c02)) + (A14*L2*b23)/l23^2 - (2*A12*L1*c02)/(b02*c23 - b23*c02), (3*L2*(3*b02*c23 - 3*b23*c02 + 2*b02*c30 - 2*b30*c02 + b23*c30 - b30*c23 + b02*c23*mu23 - b23*c02*mu23 + b23*c30*mu23 - b30*c23*mu23))/(2*(b23*c30 - b30*c23)), (A12*L2*(3*b23 + 2*b30 + b23*mu23))/(2*(b23*c30 - b30*c23)) + (A14*L2*c23)/l23^2, (A14*L2*b23)/l23^2 - (A12*L2*(3*c23 + 2*c30 + c23*mu23))/(2*(b23*c30 - b30*c23)), -(6*L3*(b02*c23 - b23*c02))/(b02*c30 - b30*c02), (2*A12*L3*b02)/(b02*c30 - b30*c02), -(2*A12*L3*c02)/(b02*c30 - b30*c02)];
D2L23 =[ (3*L1*(b02*c23 - b23*c02 + 2*b02*c30 - 2*b30*c02 + 3*b23*c30 - 3*b30*c23 - b02*c23*mu23 + b23*c02*mu23 - b23*c30*mu23 + b30*c23*mu23))/(2*(b02*c23 - b23*c02)), (A14*L1*c23)/l23^2 - (A12*L1*(2*b02 + 3*b23 - b23*mu23))/(2*(b02*c23 - b23*c02)), (A12*L1*(2*c02 + 3*c23 - c23*mu23))/(2*(b02*c23 - b23*c02)) + (A14*L1*b23)/l23^2, (3*L1*(3*b02*c23 - 3*b23*c02 + 2*b02*c30 - 2*b30*c02 + b23*c30 - b30*c23 + b02*c23*mu23 - b23*c02*mu23 + b23*c30*mu23 - b30*c23*mu23))/(2*(b23*c30 - b30*c23)) - (6*L2*(b02*c30 - b30*c02))/(b23*c30 - b30*c23), (A12*L1*(3*b23 + 2*b30 + b23*mu23))/(2*(b23*c30 - b30*c23)) + (A14*L1*c23)/l23^2 - (2*A12*L2*b30)/(b23*c30 - b30*c23), (A14*L1*b23)/l23^2 - (A12*L1*(3*c23 + 2*c30 + c23*mu23))/(2*(b23*c30 - b30*c23)) + (2*A12*L2*c30)/(b23*c30 - b30*c23), -(6*L3*(b23*c30 - b30*c23))/(b02*c30 - b30*c02), -(2*A12*L3*b30)/(b02*c30 - b30*c02), (2*A12*L3*c30)/(b02*c30 - b30*c02)];
 
% Extract from formulation (4.58), FEM - Zienkiewicz vol.2
 
 for i=1:edof
   DNL=zeros(3,3);
   DNL=[D2L11(i) D2L12(i) D2L13(i);...
       D2L12(i) D2L22(i) D2L23(i);...
       D2L13(i) D2L23(i) D2L33(i)];
%  Second derivatives of shape function coresponding to x,y in (4.58), pp.132
   DNxy=Q*DNL*Q';         
   B(1,i)=DNxy(1,1);
   B(2,i)=DNxy(2,2);
   B(3,i)=2*DNxy(1,2);
 end   
B=(1/A12^2)*B;
% Formula (27.32) in IFEM,chapter 27.
ke=ke+0.5*A12*we*B'*H*B;

end

