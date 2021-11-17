% Programming for Stiffened Plate bending using conforming HCT element
% Static Analysis of  plate
% Problem : To find the maximum bedning of plate when uniform transverse
% pressure is applied. 
% Boundary conditions are used, clamped, simply supported
%--------------------------------------------------------------------------
% Code is written by : Tran Van Nha                                       |
%                   Math & Compute Science                                |
%                                                                         |
% E-mail : tvnha108@gmail.com                                             |
%--------------------------------------------------------------------------
%|________________________________________________________________________|
%|     Edge  :3Node [w,Qx,Qy]            Model:4Node-12Dof                |
%|     Center:1Node [w]                                                   |
%|                                                                        |
%|                                         3-node (2)                     |
%|                                           o                            |
%|                                          / \                           |
%|                                                                        |
%|								L1=0      /		\   L3=0                  |
%|                                                                        |
%|										/ Sub-T_1 \                       |
%|                                                                        |
%|                        0-node  (3) o-------------o 2-node (1)          |
%|                                            L2=0                        |
%|------------------------------------------------------------------------|
%|________________________________________________________________________|
clear all
tic
format long
global sdof nel nnel nnode Lx Hy t nx ny Poisson Dmat
global gcoord ele_nods H ndof edof Isx Isy D
%----------------------------------------------%
%             1. PREPROCESSOR PHASE            %
%----------------------------------------------%
%------------------------------------------%
%  1.1 Input the initial data of example   % 
%------------------------------------------%
% nx=input('Number of element in x-axis (8,12,16,20,24,...) = ');    
% ny=input('Number of element in y-axis (8,12,16,20,24,...) = ');
nx=6;                 % Number of element in x-axis
ny=6;                 % Number of element in y-axis

fprintf('Choose example \n')
fprintf('1. Single Stiffened_ Bx.Thang \n')
fprintf('2. Cross Stiffened _ M.Barik \n')
option1=input('Option (1,2 or 3) = ');

if option1 == 1
    Example_1;
    % Boundary conditions is used, simply supported
elseif option1 == 2
    % Boundary conditions is used, simply supported
    Example_2;
% elseif option1 == 3
%     Example_2;
end

%----------------------------------------------%
%  1.2 Compute necessary data from input data  % 
%----------------------------------------------%
nel=2*nx*ny;          % Total number of elements in system
nnel=3;               % Number of nodes per element
ndof=3;               % Number of dofs per node
edof=nnel*ndof;       % Degrees of freedom per element
nnode=(nx+1)*(ny+1);  % Total number of nodes in system
sdof=nnode*3;         % Total degrees of freedom in system
D1=(Emodule*t^3)/(12*(1-Poisson^2));
Dmat=[1 Poisson 0;Poisson 1 0;0 0 (1-Poisson)/2];  % Matrix of material constants
H = D1*Dmat;                    % The Hooke matrix of Krichhoff's plate bending

[gcoord,ele_nods,bcdof,bcval]=get_initdata_HCT_SP;
%   gcoord     = matrix containing coordinate values of each node
%   ele_nods   = matrix containing nodes of each element 
%   bcdof      = vector containing dofs associated with boundary conditions
%   bcval      = vector containing boundary condition values associated with dofs in bcdof

% Showmesh_Triangle(gcoord,ele_nods)

%-----------------------------------------------%
%             2. SOLUTION PHASE                 %
%-----------------------------------------------%
%-----------------------------------------------%
%      2.1 Compute the stiffness matrix         % 
%-----------------------------------------------%

[ Kp,Fp ] = cal_Kp_Fp(Load);

%-----------------------------------------------%
%      2.2 Adding Stiffened                     %
%-----------------------------------------------%

[Kp]=Assembling_Stiffened(option1,Emodule,Kp);

%-----------------------------------------------%
%       2.2 Apply the boundary condition        %
%-----------------------------------------------%
[Kp Fp]=apply_bcdof_SP(Kp,Fp,bcdof,bcval);
%-----------------------------------------------------------------------%
% 2.3 Solve the system of equations to obtain nodal displacement vector %
%-----------------------------------------------------------------------%
% %-----------------------------------------------------------------------%
% % 2.3 Solve the system of equations to obtain nodal displacement vector %
% %-----------------------------------------------------------------------%
% %-----------------------------------------------%
% %              3. POSTPROCESSOR PHASE           %
% %-----------------------------------------------%
U=Kp\Fp;
m=1:3:sdof;
% HAU XU LY
% Do vong cua tam la w:
w=U(m);
clear m;
fh = figure ;
set(fh,'name','Preprocessing for FEA','numbertitle','off','color','w') ;
x=0:Lx/nx:Lx;
y=0:Hy/ny:Hy;
for j=1:ny+1
    for i=1:nx+1
        node=(j-1)*(nx+1)+i;
        z(i,j)=U(node*3-2);
    end
end
surf(x,y,z);
title('Finite Element Mesh of Plate') ;
% wmax=max(w)
c1 = find(gcoord(:,2)==max(gcoord(:,2))/2);  % at y = Hy/2 (along Y-axes);
c2 = find(gcoord(:,1)==max(gcoord(:,1))/2);  % at x = Lx/2 (along X-axes);
cc=intersect(c1,c2);% taam
wc_fem=U(cc*3-2)
wcss=wc_fem*100*D1/(Load*(Lx)^4)
toc