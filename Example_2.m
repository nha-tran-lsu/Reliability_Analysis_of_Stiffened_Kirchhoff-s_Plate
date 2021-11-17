% Barick
Poisson=0.3;          % Poisson's ratio
Emodule=2.0684*10^5;  % Young elastic modulus  
Load=6.89*10^-2;      % Uniform load
% Plate
Lx=762;               % Length of the Plate along X-axes
Hy=1524;              % High of the Plate along Y-axes
t=6.35;               % Thickness of plate
% Beam // ox
bx=12.7;    tx=127;     ISy= bx*tx^3/12;    %ex=tx/2;    Ax=bx*tx;   ISy= bx*tx^3/12;%+ex^2*Ax;
% Beam // oy
by=12.7;    ty=76.2;    ISx=by*ty^3/12;     %ey=ty/2;    Ay=by*ty;   ISx=by*ty^3/12;%+ey^2*Ay;
