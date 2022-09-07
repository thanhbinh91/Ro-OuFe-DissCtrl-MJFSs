function [nx, nu, nw, ny, nz, s, r, A, B, E, C, D, G, H, Jj, Pi] = SysParas
%% Dimension 
nx = 2;
nu = 1; 
nw = 1; 
ny = 1; 
nz = 1; 
s  = 3; 
r  = 2; 

%% MJFS
A  = zeros(nx,nx,s,r);
B  = zeros(nx,nu,s,r);
E  = zeros(nx,nw,s,r);
C  = zeros(ny,nx,s,r);
D  = zeros(ny,nw,s,r);
G  = zeros(nz,nx,s,r);
H  = zeros(nz,nu,s,r);
Jj = zeros(nz,nw,s,r);

%% Markov chain desription
Ns = [1 2 3];


Pi  = [0.8       0.1    0.1;
       0.2       0.7    0.1;
       0.5       0.2    0.3];

Ts = 0.5;

M = [0.75 1.5 2];
J = [1    2   2.5];

cv = 2;
l  = 0.5;
ga = 9.81;
eta = 0.01/pi;


%% Robot arm in discrete-time Markovian Fuzzy model


for g = Ns
    A(:,:,g,1) = [ 1                      Ts;
                  -Ts*M(g)*ga*l/J(g)      1-Ts*cv/J(g)];
              
    A(:,:,g,2) = [ 1                      Ts;
                  -eta*Ts*M(g)*ga*l/J(g)  1-Ts*cv/J(g)];
              
    B(:,:,g,1) = [0; Ts/J(g)];
    B(:,:,g,2) = B(:,:,g,1);
    
    E(:,:,g,1) = [0; Ts];
    E(:,:,g,2) = E(:,:,g,1);
    
    C(:,:,g,1) = [1 0];
    C(:,:,g,2) = C(:,:,g,1);
    
    D(:,:,g,1) = 0.1;
    D(:,:,g,2) = D(:,:,g,1);
    
    G(:,:,g,1) = [1 0];
    G(:,:,g,2) = G(:,:,g,1);
    
    H(:,:,g,1) = 0.1;
    H(:,:,g,2) = H(:,:,g,1);
end
end