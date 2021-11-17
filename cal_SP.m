function [ wmax ] = cal_SP(Emodule,Load )
% Function for Krichhoff's plate bending using conforming HCT element
% Static Analysis of  plate
% Problem : To find the maximum bedning of plate when uniform transverse
% pressure is applied. 
% Boundary conditions is used, clamped
%--------------------------------------------------------------------------
% Code is written by : Tran Van Nha                                          |
%                   Math & Compute Science                                |
%                                                                         |
% E-mail : tvnha108@gmail.com                                             |
%--------------------------------------------------------------------------
%|________________________________________________________________________|
%|     Edge  :3Node [w,Qx,Qy]            Model:4Node-12Dof                |
%|     Center:1Node [w]                                                   |
%|                                                                        |
%|________________________________________________________________________|
%--------------------------------------------------------------------------
% clc
% clear all
% tic
format long
global sdof nel nnel nnode Lx Hy t nx ny Poisson D Dmat
global gcoord ele_nods H ndof edof ISy ISx %Emodule Load
%----------------------------------------------%
%             1. PREPROCESSOR PHASE            %
%----------------------------------------------%
%------------------------------------------%
%  1.1 Input the initial data of example   % 
%------------------------------------------%
% nx=input('Number of element in x-axis (8,12,16,20,24,...) = ');    
% ny=input('Number of element in y-axis (8,12,16,20,24,...) = ');
nx=16;                 % Number of element in x-axis
ny=16;                 % Number of element in y-axis

Poisson=0.3;          % Poisson's ratio
% Emodule=1092000;      % Young elastic modulus  
% Load=1;               % Uniform load
Lx=1;                 % Length of the Plate along X-axes
Hy=1;                 % High of the Plate along Y-axes
t=0.001;              % Thickness of plate
ISy=10^-12;
ISx=10^-12;
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
[ Ksx ] = cal_Ksx(Emodule);
[ Ksy ] = cal_Ksy(Emodule);
% Assembling Ksx
[ Tsx, Tsy ] = Assembling_Ksx_Ksy;
Kp = Kp + Tsx'*Ksx*Tsx + Tsy'*Ksy*Tsy ;
% Kp=Kp+Tsx'*Ksx*Tsx;
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
% Chuan hoa Wc de so sanh ket qua voi tai lieu thay Trung
wmax=max(w);
wcss=wmax*100*D1/(Load*(Lx)^4);
wcssa=0.1267;
% toc
end

