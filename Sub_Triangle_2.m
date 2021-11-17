clear all
clc
syms L1 L2 L3 b10 c10 b03 c03 b31 c31 l10 l03 l31 mu10 mu03 mu31 A22 A24
%J=[1 x1 y1;1 x2 y2; 1 x3 y3];
% A22=2*area_Sub_triangle_2
% o: The centroid of triangle element

%                                          3-node (1)
%                                            o
%                                           / \

%								L3=0      /		\   L2=0

%										/ Sub-T_2 \

%                         1-node  (2) o-------------o 0-node (3)

% x0=(x1+x2+x3)/3;   y0=(y1+y2+y3)/3;
% b10= y1-y0;                 b03=y0-y3;                  b31=y3-y1;
% c10= x1-x0;                 c03=x0-x3;                  c31=x3-x1;
% l10=sqrt(b10^2+c10^2);      l03=sqrt(b03^2+c03^2);      l31=sqrt(b31^2+c31^2);
% mu10=(l31^2-l03^2)/l10^2;   mu03=(l10^2-l31^2)/l03^2;   mu31=(l03^2-l10^2)/l31^2;
% A22=b10*(-c03)-(-c10)*b03;      % A22=2*area_Sub_triangle_2
% A24=2*A22
C2=-c31/l31;        % Cos(phi_2)
S2=-b31/l31;        % Sin(phi_2)

P=[L1^3 L2^3 L3^3 L1^2*L2 L1^2*L3 L2^2*L3 L2^2*L1 L3^2*L1 L3^2*L2 L1*L2*L3];

S =1/A22*[...
    [ A22,    0,    0,    0];...
    [   0,  c10,  c03,  c31];...
    [   0,  b10,  b03,  b31]];

Q=[S            zeros(3,4)  zeros(3,4)  zeros(3,1);...
   zeros(3,4)   S           zeros(3,4)  zeros(3,1);...
   zeros(3,4)   zeros(3,4)  S           zeros(3,1);...
   zeros(1,4)   zeros(1,4)  zeros(1,4)  1         ];

M=[ subs(P,{L1,L2,L3},{1,0,0});...
    subs(diff(P,L1),{L1,L2,L3},{1,0,0});...
    subs(diff(P,L2),{L1,L2,L3},{1,0,0});...
    subs(diff(P,L3),{L1,L2,L3},{1,0,0});...
    subs(P,{L1,L2,L3},{0,1,0});...
    subs(diff(P,L1),{L1,L2,L3},{0,1,0});...
    subs(diff(P,L2),{L1,L2,L3},{0,1,0});...
    subs(diff(P,L3),{L1,L2,L3},{0,1,0});...
    subs(P,{L1,L2,L3},{0,0,1});...
    subs(diff(P,L1),{L1,L2,L3},{0,0,1});...
    subs(diff(P,L2),{L1,L2,L3},{0,0,1});...
    subs(diff(P,L3),{L1,L2,L3},{0,0,1});...
    l31/A24*(subs(diff(P,L1)+diff(P,L2)-2*diff(P,L3)+mu31*(diff(P,L2)-diff(P,L1)),{L1,L2,L3},{1/2,1/2,0}))];
C=Q*M;
N=P*inv(C);
N1= N(1);    N2= N(2) + (C2*N(10))/2;    N3= N(3) + (N(10)*S2)/2;
N4= N(4);    N5= N(5) + (C2*N(10))/2;    N6= N(6) + (N(10)*S2)/2;
N7= N(7);    N8= N(8);                   N9= N(9);
SN2=[N1 N2 N3 N4 N5 N6 N7 N8 N9];
D2L11=diff(diff(SN2,L1),L1)
D2L22=diff(diff(SN2,L2),L2)
D2L33=diff(diff(SN2,L3),L3)
D2L12=diff(diff(SN2,L1),L2)
D2L13=diff(diff(SN2,L1),L3)
D2L23=diff(diff(SN2,L2),L3)
%w_x =-dw/dy
NN2=subs(SN2,L3,1-L1-L2);
fe2=int(int(NN2,L2,0,1-L1),L1,0,1)

SN2e=[N4 N5 N6 0 0 0 N1 N2 N3];
SN2o=[N7 N8 N9];

% B2e9=d(SN2e)/dn|9
B2e9 = l03/A24*(subs(diff(SN2e,L3)+diff(SN2e,L1)-2*diff(SN2e,L2)+mu03*(diff(SN2e,L1)-diff(SN2e,L3)),{L1,L2,L3},{1/2,0,1/2}))
% B2o9=d(SN2o)/dn|9
B2o9 = l03/A24*(subs(diff(SN2o,L3)+diff(SN2o,L1)-2*diff(SN2o,L2)+mu03*(diff(SN2o,L1)-diff(SN2o,L3)),{L1,L2,L3},{1/2,0,1/2}))
% B2e7=d(SN2e)/dn|7
B2e7 = l10/A24*(subs(diff(SN2e,L2)+diff(SN2e,L3)-2*diff(SN2e,L1)+mu10*(diff(SN2e,L3)-diff(SN2e,L2)),{L1,L2,L3},{0,1/2,1/2}))
% B2o7=d(SN2o)/dn|7
B2o7 = l10/A24*(subs(diff(SN2o,L2)+diff(SN2o,L3)-2*diff(SN2o,L1)+mu10*(diff(SN2o,L3)-diff(SN2o,L2)),{L1,L2,L3},{0,1/2,1/2}))
% ///////////////////////////////////////////////////////////////////////////////////////