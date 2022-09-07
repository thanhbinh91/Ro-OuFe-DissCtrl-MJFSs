function [Xg, Xg_, Pg, Fg, Lg, tmin] = Theorem3_S1(mode, al, beta)

[nx, nu, nw, ny, nz, s, r, A, B, E, C, D, G, H, J, Pi] = SysParas;

%% Variables Declaration
nlmis = 1;
setlmis([]);
disp('=== Theorem 3: Step 1 LMIs solvers ===')
if mode == 1
    Q = -eye(nz); Q1 = chol(-Q); nq = nz;
    S =  zeros(nw,nz);
%   H infinity performance: Minimize beta^2');
%   R = beta^2 + beta,'); 
elseif mode == 2
    Q = 0*eye(nz); Q1 = 0*eye(nz); nq = nz;
    S = eye(nw,nz);
%   Passivity performance: Minimize beta
else
    Q = -0.01*eye(nz); Q1 = chol(-Q); nq = nz;
    S =  0.2*ones(nw,nz);
    R =  5*eye(nw);
%   Dissipative performance: Q, S, R and Minimize -beta
end




L   = zeros(s,r);
P   = zeros(s,r,r);
F   = zeros(s,r);
W   = zeros(s,r);
U   = zeros(s,r,r);

X  = lmivar(1, [2*nx  1]);
X_ = lmivar(1, [2*nx  1]);

for g = 1:s
    for i = 1:r
    for j = 1:r
        P(g,i,j) = lmivar(1, [2*nx  1]);
    end
    end
end


for g = 1:s
for i = 1:r
    F(g,i) = lmivar(2, [nu  nx]);
    L(g,i) = lmivar(2, [nx  ny]);
    W(g,i) = lmivar(2, [2*nx+nw 2*nx+nw+nz]);
    for j = 1:r
        U(g,i,j) = lmivar(1, [2*nx+nw 1]);
    end
end
end


%% Lmi condition
for g = 1:s
for i = 1:r
    lmiterm([-nlmis  1     1     P(g,i,i)],  1,        1);
    nlmis = nlmis + 1;
    %
    lmiterm([-nlmis  1     1     X],      1,        1);
    for h = 1:s
    lmiterm([-nlmis  1     1     P(h,i,i)], -Pi(g,h),        1);
    end
    nlmis = nlmis + 1;
    %
    Xi_(1,g,i,i);
    nlmis = nlmis + 1;
for j = 1:r
    if j ~= i
    lmiterm([-nlmis  1     1     P(g,i,i)],  1/(r-1),  1);
    lmiterm([-nlmis  1     1     P(g,i,j)],  1,        1/2); 
    lmiterm([-nlmis  1     1     P(g,j,i)],  1,        1/2); 
    nlmis = nlmis + 1;
    %
    lmiterm([-nlmis  1     1     X],      r/(r-1),  1);
    for h = 1:s
    lmiterm([-nlmis  1     1     P(h,i,i)], -Pi(g,h)/(r-1),  1);
    lmiterm([-nlmis  1     1     P(h,i,j)], -Pi(g,h),        1/2); 
    lmiterm([-nlmis  1     1     P(h,j,i)], -Pi(g,h),        1/2); 
    end
    nlmis = nlmis + 1;        
    %
    Xi_(2/(r-1),g,i,i);
    Xi_(1,g,i,j);
    Xi_(1,g,j,i);
    nlmis = nlmis + 1;
    end
end
end
end


%%%%%%%%%%%%%%% X X_ = I %%%%%%%%%%%%%%%%%%%%

lmiterm([-nlmis  1     1     X],  1,        1);
lmiterm([-nlmis  2     2     X_], 1,        1);
lmiterm([-nlmis  2     1     0],     eye(2*nx));
nlmis = nlmis + 1;



lmis = getlmis;
% fprintf('number of matrix variables = %d\n',matnbr(lmis));
% fprintf('number of LMIs = %d\n', lminbr(lmis));
% fprintf('number of variables = %d\n', decnbr(lmis));
options = [1e-3,200,0,0,1];



[tmin,xopt] = feasp(lmis,options);

Pg  = zeros(2*nx,2*nx,g,r,r);
Fg  = zeros(nu,nx,g,r);
Lg  = zeros(nx,ny,g,r);

Xg  = dec2mat(lmis, xopt, X);
Xg_ = dec2mat(lmis, xopt, X_);

if tmin < 0
    fprintf('Step 1: feasible\n\n');
    for g = 1:s    
    for i = 1:r
    Fg(:,:,g,i)  = dec2mat(lmis, xopt, F(g,i));
    Lg(:,:,g,i)  = dec2mat(lmis, xopt, L(g,i));
        for j = 1:r
            Pg(:,:,g,i,j)  = dec2mat(lmis, xopt, P(g,i,j));
        end
    end
    end

else
    disp('Step 1: Infeasible')
end




%% ================Nested Functions================== %%
function  Xi_(co,g,i,j)
%     nN = max(size(Ng{g}));
    
    Om1  = [eye(nx)      zeros(nx)    zeros(nx,nw);
            zeros(nx)    eye(nx)      zeros(nx,nw);
            zeros(nw,nx) zeros(nw,nx) eye(nw);
            zeros(nq,nx) zeros(nq,nx) zeros(nq,nw);
            zeros(2*nx)  zeros(2*nx,nw)];
    Om1  =  Om1';
    
    Om2  = [zeros(nw,nx) zeros(nw,nx) eye(nw)      zeros(nw,nq); 
            zeros(nq,nx) zeros(nq,nx) zeros(nq,nw) eye(nq); 
            zeros(nx,nx+nx+nw+nq);
            zeros(nx,nx+nx+nw+nq)];
    Omm2 = [zeros(nw,nx) zeros(nw,nx);
            zeros(nq,nx) zeros(nq,nx);
            eye(nx)      zeros(nx,nx);
            zeros(nx,nx) eye(nx)];
    Om2  = [Om2 Omm2];

    e = cell(4,1);
    e{1} = [eye(2*nx),          zeros(2*nx, nw+nq+2*nx)];
    e{2} = [zeros(nw,2*nx),     eye(nw), zeros(nw,nq+2*nx)];
    e{3} = [zeros(nq,2*nx+nw),  eye(nq), zeros(nq,2*nx)];
    e{4} = [zeros(2*nx,2*nx+nw+nq), eye(2*nx)];

    
    %% (1,1) - Xi0
    if mode == 1
    lmiterm([nlmis  1  1  0],          -co*beta^2*e{2}'*e{2});
    elseif mode == 2
    lmiterm([nlmis  1  1  0],          -co*beta*e{2}'*e{2});
    else
    lmiterm([nlmis  1  1  0],           co*beta*e{2}'*e{2});
    lmiterm([nlmis  1  1  0],          -co*e{2}'*R*e{2});
    end
    %
    lmiterm([nlmis  1  1  0],          -co*e{3}'*e{3});
    %
    lmiterm([nlmis  1  1  X_],      -co*e{4}',   e{4});
    %
    %================= Up1 ====================
    Xi1 = [-G(:,:,g,i)'*S',  G(:,:,g,i)'*Q1',  zeros(nx),     A(:,:,g,i)';
           -G(:,:,g,i)'*S',  G(:,:,g,i)'*Q1',  zeros(nx),     A(:,:,g,i)';
           -J(:,:,g,i)'*S',  J(:,:,g,i)'*Q1',  zeros(nw,nx),  E(:,:,g,i)'];
       
    lmiterm([nlmis  1  1  0],           co*(Om1'*Xi1*Om2 + Om2'*Xi1'*Om1));
    %
    %================= Up2 ====================
    p1 = [eye(2*nx),        zeros(2*nx,nw),   zeros(2*nx,nq), zeros(2*nx,2*nx)];
    p2 = zeros(nx,2*nx+nw+nq);  p2 = [p2,  eye(nx),     zeros(nx)];
    p3 = zeros(nx,2*nx+nw+nq);  p3 = [p3,  zeros(nx),   eye(nx)];
    p4 = [eye(nx), zeros(nx,nx+nw+nq+2*nx)];
    
    
%     lmiterm([nlmis  1  1  P(g,i)],    -co*p1',               p1);
    lmiterm([nlmis  1  1  0],          co*(p2'*A(:,:,g,i)*p4 + p4'*A(:,:,g,i)'*p2));
    lmiterm([nlmis  1  1  0],         -co*(p3'*A(:,:,g,i)*p4 + p4'*A(:,:,g,i)'*p3));
    
    %================= Up3 ====================
    c1 = [eye(nw)      zeros(nw,nz) zeros(nw,nx) zeros(nw,nx)];
    c2 = [zeros(nz,nw) eye(nz)      zeros(nz,nx) zeros(nz,nx)];
    c3 = [zeros(nx,nw) zeros(nx,nz) eye(nx)      zeros(nx)];
    c4 = [zeros(nx,nw) zeros(nx,nz) zeros(nx)    eye(nx)];
    l1 = [eye(nx)      zeros(nx)    zeros(nx,nw)];
    l2 = [zeros(nx)    eye(nx)      zeros(nx,nw)];
    l3 = [zeros(nw,nx) zeros(nw,nx) eye(nw)];
    
    lmiterm([nlmis  1  1 -F(g,j)],     -co*Om1'*l1',             H(:,:,g,i)'*S'*c1*Om2, 's');
    lmiterm([nlmis  1  1 -F(g,j)],      co*Om1'*l1',             H(:,:,g,i)'*Q1'*c2*Om2,'s');
    lmiterm([nlmis  1  1 -L(g,j)],      co*Om1'*l1'*C(:,:,g,i)', c3*Om2, 's');
    lmiterm([nlmis  1  1 -F(g,j)],      co*Om1'*l1',             B(:,:,g,i)'*c4*Om2, 's');
    lmiterm([nlmis  1  1 -L(g,j)],     -co*Om1'*l1'*C(:,:,g,i)', c4*Om2, 's');%%%%
    lmiterm([nlmis  1  1 -L(g,j)],      co*Om1'*l2'*C(:,:,g,i)', c3*Om2, 's');
    lmiterm([nlmis  1  1 -L(g,j)],     -co*Om1'*l2'*C(:,:,g,i)', c4*Om2, 's');
    lmiterm([nlmis  1  1 -L(g,j)],      co*Om1'*l3'*D(:,:,g,i)', c3*Om2, 's');
    lmiterm([nlmis  1  1 -L(g,j)],     -co*Om1'*l3'*D(:,:,g,i)', c4*Om2, 's');%%%%

    %================= Up4 ====================
    lmiterm([nlmis  1  1  P(g,i,j)],   -co*p1',                 p1);
    lmiterm([nlmis  1  1  F(g,j)],      co*p2'*B(:,:,g,i),      p4, 's');
    lmiterm([nlmis  1  1  L(g,j)],     -co*p2',                 C(:,:,g,i)*p4, 's');
    lmiterm([nlmis  1  1  F(g,j)],     -co*p3'*B(:,:,g,i),      p4, 's');
    lmiterm([nlmis  1  1  L(g,j)],      co*p3',                 C(:,:,g,i)*p4, 's');
    %
    for kk = 1:r
    lmiterm([nlmis  1  1  U(g,kk,i)],  -co*al^2*Om1', Om1);    
    end
        

    %% (2,1)
    cR = 1; cC = 1;
    for l = 1:r
    cR = cR + 1; cC = cC + 1;
    Xi1 = [-G(:,:,g,l)'*S',  G(:,:,g,l)'*Q1',  zeros(nx),    A(:,:,g,l)';
           -G(:,:,g,l)'*S',  G(:,:,g,l)'*Q1',  zeros(nx),    A(:,:,g,l)';
           -J(:,:,g,l)'*S',  J(:,:,g,l)'*Q1',  zeros(nw,nx), E(:,:,g,l)'];      
    lmiterm([nlmis  cR  1  0],        co*Xi1*Om2);    
    
    lmiterm([nlmis  cR  1 -F(g,i)],    -co*l1',                H(:,:,g,l)'*S'*c1*Om2);
    lmiterm([nlmis  cR  1 -F(g,i)],     co*l1',                H(:,:,g,l)'*Q1'*c2*Om2);
    lmiterm([nlmis  cR  1 -L(g,i)],     co*l1'*C(:,:,g,l)',    c3*Om2);
    lmiterm([nlmis  cR  1 -F(g,i)],     co*l1',                B(:,:,g,l)'*c4*Om2);
    lmiterm([nlmis  cR  1 -L(g,i)],    -co*l1'*C(:,:,g,l)',    c4*Om2);%%%%
    lmiterm([nlmis  cR  1 -L(g,i)],     co*l2'*C(:,:,g,l)',    c3*Om2);
    lmiterm([nlmis  cR  1 -L(g,i)],    -co*l2'*C(:,:,g,l)',    c4*Om2);
    lmiterm([nlmis  cR  1 -L(g,i)],     co*l3'*D(:,:,g,l)',    c3*Om2);
    lmiterm([nlmis  cR  1 -L(g,i)],    -co*l3'*D(:,:,g,l)',    c4*Om2);%%%%    
    
    lmiterm([nlmis  cR  1  W(g,i)],     co,                      Om2);
    %% (2,2)
    lmiterm([nlmis  cR  cC U(g,l,i)],   co,                      1);
    end
end
        
end