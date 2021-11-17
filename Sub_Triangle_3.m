clear all
clc
syms L1 L2 L3 b20 c20 l20 mu20 b01 c01 l01 mu01 b12 c12 l12 mu12 A32 A34
%J=[1 x1 y1;1 x2 y2; 1 x3 y3];
% A32=2*area_Sub_triangle_3
% o: The centroid of triangle element

%                                          o-node (3)
%                                            o
%                                           / \

%								L2=0      /		\   L1=0

%										/ Sub-T_1 \

%                         1-node  (1) o-------------o 3-node (2)
%                                            L3=0 
% x0=(x1+x2+x3)/3;   y0=(y1+y2+y3)/3;
% b20= y2-y0;                 b01=y0-y1;                  b12=y1-y2;
% c20= x2-x0;                 c01=x0-x1;                  c12=x1-x2;
% l20=sqrt(b20^2+c20^2);      l01=sqrt(b01^2+c01^2);      l12=sqrt(b12^2+c12^2);
% mu20=(l12^2-l01^2)/l20^2;   mu01=(l20^2-l12^2)/l01^2;   mu12=(l01^2-l20^2)/l12^2;
% A32=b20*(-c01)-(-c20)*b01;      % A22=2*area_Sub_triangle_3
% A34=2*A32
C3=-c12/l12;        % Cos(phi_3)
S3=-b12/l12;        % Sin(phi_3)

P=[L1^3 L2^3 L3^3 L1^2*L2 L1^2*L3 L2^2*L3 L2^2*L1 L3^2*L1 L3^2*L2 L1*L2*L3];

S =1/A32*[...
    [ A32,    0,    0,    0];...
    [   0,  c20,  c01,  c12];...
    [   0,  b20,  b01,  b12]];

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
    l12/A34*(subs(diff(P,L1)+diff(P,L2)-2*diff(P,L3)+mu12*(diff(P,L2)-diff(P,L1)),{L1,L2,L3},{1/2,1/2,0}))];
C=Q*M;
N=P*inv(C);
N1= N(1);    N2= N(2) + (C3*N(10))/2;    N3= N(3) + (N(10)*S3)/2;
N4= N(4);    N5= N(5) + (C3*N(10))/2;    N6= N(6) + (N(10)*S3)/2;
N7= N(7);    N8= N(8);                   N9= N(9);
SN3=[N1 N2 N3 N4 N5 N6 N7 N8 N9];
D2L11=diff(diff(SN3,L1),L1);
D2L22=diff(diff(SN3,L2),L2);
D2L33=diff(diff(SN3,L3),L3);
D2L12=diff(diff(SN3,L1),L2);
D2L13=diff(diff(SN3,L1),L3);
D2L23=diff(diff(SN3,L2),L3);
%w_x =-dw/dy
NN3=subs(SN3,L3,1-L1-L2);
fe3=int(int(NN3,L2,0,1-L1),L1,0,1)

SN3e=[N1 N2 N3 N4 N5 N6 0 0 0];
SN3o=[N7 N8 N9];
% B3e7=d(SN3e)/dn|7
B3e7 = l01/A34*(subs(diff(SN3e,L3)+diff(SN3e,L1)-2*diff(SN3e,L2)+mu01*(diff(SN3e,L1)-diff(SN3e,L3)),{L1,L2,L3},{1/2,0,1/2}))
% B3o7=d(SN3o)/dn|7
B3o7 = l01/A34*(subs(diff(SN3o,L3)+diff(SN3o,L1)-2*diff(SN3o,L2)+mu01*(diff(SN3o,L1)-diff(SN3o,L3)),{L1,L2,L3},{1/2,0,1/2}))
% B3e8=d(SN3e)/dn|8
B3e8 = l20/A34*(subs(diff(SN3e,L2)+diff(SN3e,L3)-2*diff(SN3e,L1)+mu20*(diff(SN3e,L3)-diff(SN3e,L2)),{L1,L2,L3},{0,1/2,1/2}))
% B3o8=d(SN3o)/dn|8
B3o8 = l20/A34*(subs(diff(SN3o,L2)+diff(SN3o,L3)-2*diff(SN3o,L1)+mu20*(diff(SN3o,L3)-diff(SN3o,L2)),{L1,L2,L3},{0,1/2,1/2}))
% ///////////////////////////////////////////////////////////////////////////////////////