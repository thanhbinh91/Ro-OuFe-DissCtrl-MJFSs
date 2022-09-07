function [nx, nu, nw, ny, nz, r, A, B, E, C, D, G, H, J] = SysParas
%% Dimension 
nx = 3;
nu = 1; 
nw = 1; 
ny = 3; 
nz = 1; 
r  = 2;

%% MJFS
A  = zeros(nx,nx,r);
B  = zeros(nx,nu,r);
E  = zeros(nx,nw,r);
C  = zeros(ny,nx,r);
D  = zeros(ny,nw,r);
G  = zeros(nz,nx,r);
H  = zeros(nz,nu,r);
J  = zeros(nz,nw,r);




%% Markov chain desription
%% 
    l = 2.8;
    L = 5.5;
    v = -1;
    Ts  =  2;
    g = 0.01/pi;

    A(:,:,1) = [1-v*Ts/L    0      0;
                v*Ts/L      1      0;
               (v*Ts)^2/L   v*Ts   1];
              
    A(:,:,2) = [1-v*Ts/L       0        0;
                v*Ts/L         1        0;
                g*(v*Ts)^2/L   g*v*Ts   1];
              
    B(:,:,1) = [v*Ts/l; 0; 0];
    B(:,:,2) = B(:,:,1);
    
    E(:,:,1) = [0; 0.2; 0.1];
    E(:,:,2) = E(:,:,1);
    
    C(:,:,1) = [1 0 1;
                0 2 1;
                1 2 2];
            
    C(:,:,2) = [1 0 1;
                0 1 1;
                1 1 1];
    
    D(:,:,1) = 0;
    D(:,:,2) = D(:,:,1);
    
    G(:,:,1)  = [ 0.1 0 0];
    G(:,:,2)  = [-0.1 0 0];
    
    H(:,:,1)  = -0.1;
    H(:,:,2)  = -0.1;
     
    J(:,:,1) = 3;
    J(:,:,2) = -3;
   

end