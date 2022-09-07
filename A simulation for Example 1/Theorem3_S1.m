function [Xg, Xg_, Pg, Fg, Lg, tmin] = Theorem3_S1(mode, al, beta)
[nx, nu, nw, ny, nz, r, A, B, E, C, D, G, H, J] = SysParas;

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
    Q = -0.36*eye(nz); Q1 = chol(-Q); nq = nz;
    S = -2*ones(nw,nz);
    R =  2*eye(nw);
%   Dissipative performance: Q, S, R and Minimize -beta
end




L   = zeros(r);
P   = zeros(r,r);
F   = zeros(r);
W   = zeros(r);
U   = zeros(r,r);

X  = lmivar(1, [2*nx  1]);
X_ = lmivar(1, [2*nx  1]);

for i = 1:r
    for j = 1:r
        P(i,j) = lmivar(1, [2*nx  1]);
    end
end



for i = 1:r
    F(i) = lmivar(2, [nu  nx]);
    L(i) = lmivar(2, [nx  ny]);
    W(i) = lmivar(2, [2*nx+nw 2*nx+nw+nz]);
    for j = 1:r
        U(i,j) = lmivar(1, [2*nx+nw 1]);
    end
end



%% Lmi condition
for i = 1:r
    lmiterm([-nlmis  1     1     P(i,i)],  1,        1); 
    nlmis = nlmis + 1;
    %
    lmiterm([-nlmis  1     1     X],       1,        1);
    lmiterm([-nlmis  1     1     P(i,i)], -1,        1);
    nlmis = nlmis + 1;
    %
    Xi_(1,i,i);
    nlmis = nlmis + 1;
for j = 1:r
    if j ~= i
    lmiterm([-nlmis  1     1     P(i,i)],  1/(r-1),  1);
    lmiterm([-nlmis  1     1     P(i,j)],  1,        1/2); 
    lmiterm([-nlmis  1     1     P(j,i)],  1,        1/2); 
    nlmis = nlmis + 1;
    %
    lmiterm([-nlmis  1     1     X],       r/(r-1),  1);
    lmiterm([-nlmis  1     1     P(i,i)], -1/(r-1),  1);
    lmiterm([-nlmis  1     1     P(i,j)], -1,        1/2); 
    lmiterm([-nlmis  1     1     P(j,i)], -1,        1/2); 
    nlmis = nlmis + 1;        
    %
    Xi_(2/(r-1),i,i);
    Xi_(1,i,j);
    Xi_(1,j,i);
    nlmis = nlmis + 1;
    end
end
end


%%%%%%%%%%%%%%% Xg Xg_ = I %%%%%%%%%%%%%%%%%%%%
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

Pg  = zeros(2*nx,2*nx,r,r);
Xg  = zeros(2*nx,2*nx);
Xg_ = zeros(2*nx,2*nx);
Fg  = zeros(nu,nx,r);
Lg  = zeros(nx,ny,r);

if tmin < 0
    fprintf('Step 1: feasible\n\n');
    
    Xg(:,:)  = dec2mat(lmis, xopt, X);
    Xg_(:,:) = dec2mat(lmis, xopt, X_);
    for i = 1:r
        Fg(:,:,i)  = dec2mat(lmis, xopt, F(i));
        Lg(:,:,i)  = dec2mat(lmis, xopt, L(i));
        for j = 1:r
            Pg(:,:,i,j)  = dec2mat(lmis, xopt, P(i,j));
        end
    end
else
    disp('Step 1: Infeasible')
end




%% ================Nested Functions================== %%
function  Xi_(co,i,j)
    
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
    lmiterm([nlmis  1  1  X_],         -co*e{4}',   e{4});
    %
    %================= Up1 ====================
    Xi1 = [-G(:,:,i)'*S',  G(:,:,i)'*Q1',  zeros(nx),     A(:,:,i)';
           -G(:,:,i)'*S',  G(:,:,i)'*Q1',  zeros(nx),     A(:,:,i)';
           -J(:,:,i)'*S',  J(:,:,i)'*Q1',  zeros(nw,nx),  E(:,:,i)'];
       
    lmiterm([nlmis  1  1  0],           co*(Om1'*Xi1*Om2 + Om2'*Xi1'*Om1));
    %
    %================= Up2 ====================
    p1 = [eye(2*nx),    zeros(2*nx,nw),   zeros(2*nx,nq), zeros(2*nx,2*nx)];
    p2 = zeros(nx,2*nx+nw+nq);  p2 = [p2,  eye(nx),     zeros(nx)];
    p3 = zeros(nx,2*nx+nw+nq);  p3 = [p3,  zeros(nx),   eye(nx)];
    p4 = [eye(nx), zeros(nx,nx+nw+nq+2*nx)];
    
    
%     lmiterm([nlmis  1  1  P(i)],    -co*p1',               p1);
    lmiterm([nlmis  1  1  0],          co*(p2'*A(:,:,i)*p4 + p4'*A(:,:,i)'*p2));
    lmiterm([nlmis  1  1  0],         -co*(p3'*A(:,:,i)*p4 + p4'*A(:,:,i)'*p3));
    
    %================= Up3 ====================
    c1 = [eye(nw)      zeros(nw,nz) zeros(nw,nx) zeros(nw,nx)];
    c2 = [zeros(nz,nw) eye(nz)      zeros(nz,nx) zeros(nz,nx)];
    c3 = [zeros(nx,nw) zeros(nx,nz) eye(nx)      zeros(nx)];
    c4 = [zeros(nx,nw) zeros(nx,nz) zeros(nx)    eye(nx)];
    l1 = [eye(nx)      zeros(nx)    zeros(nx,nw)];
    l2 = [zeros(nx)    eye(nx)      zeros(nx,nw)];
    l3 = [zeros(nw,nx) zeros(nw,nx) eye(nw)];
    
    lmiterm([nlmis  1  1 -F(j)],     -co*Om1'*l1',             H(:,:,i)'*S'*c1*Om2, 's');
    lmiterm([nlmis  1  1 -F(j)],      co*Om1'*l1',             H(:,:,i)'*Q1'*c2*Om2,'s');
    lmiterm([nlmis  1  1 -L(j)],      co*Om1'*l1'*C(:,:,i)', c3*Om2, 's');
    lmiterm([nlmis  1  1 -F(j)],      co*Om1'*l1',             B(:,:,i)'*c4*Om2, 's');
    lmiterm([nlmis  1  1 -L(j)],     -co*Om1'*l1'*C(:,:,i)', c4*Om2, 's');%%%%
    lmiterm([nlmis  1  1 -L(j)],      co*Om1'*l2'*C(:,:,i)', c3*Om2, 's');
    lmiterm([nlmis  1  1 -L(j)],     -co*Om1'*l2'*C(:,:,i)', c4*Om2, 's');
    lmiterm([nlmis  1  1 -L(j)],      co*Om1'*l3'*D(:,:,i)', c3*Om2, 's');
    lmiterm([nlmis  1  1 -L(j)],     -co*Om1'*l3'*D(:,:,i)', c4*Om2, 's');%%%%

    %================= Up4 ====================
    lmiterm([nlmis  1  1  P(i,j)],   -co*p1',               p1);
    lmiterm([nlmis  1  1  F(j)],      co*p2'*B(:,:,i),      p4, 's');
    lmiterm([nlmis  1  1  L(j)],     -co*p2',               C(:,:,i)*p4, 's');
    lmiterm([nlmis  1  1  F(j)],     -co*p3'*B(:,:,i),      p4, 's');
    lmiterm([nlmis  1  1  L(j)],      co*p3',               C(:,:,i)*p4, 's');
    %
    for kk = 1:r
    lmiterm([nlmis  1  1  U(kk,i)],  -co*al^2*Om1', Om1);    
    end
        

    %% (2,1)
    cR = 1; cC = 1;
    for l = 1:r
    cR = cR + 1; cC = cC + 1;

    Xi1 = [-G(:,:,l)'*S',  G(:,:,l)'*Q1',  zeros(nx),    A(:,:,l)';
           -G(:,:,l)'*S',  G(:,:,l)'*Q1',  zeros(nx),    A(:,:,l)';
           -J(:,:,l)'*S',  J(:,:,l)'*Q1',  zeros(nw,nx), E(:,:,l)'];      
    lmiterm([nlmis  cR  1  0],        co*Xi1*Om2);    
    
    lmiterm([nlmis  cR  1 -F(i)],    -co*l1',                H(:,:,l)'*S'*c1*Om2);
    lmiterm([nlmis  cR  1 -F(i)],     co*l1',                H(:,:,l)'*Q1'*c2*Om2);
    lmiterm([nlmis  cR  1 -L(i)],     co*l1'*C(:,:,l)',    c3*Om2);
    lmiterm([nlmis  cR  1 -F(i)],     co*l1',                B(:,:,l)'*c4*Om2);
    lmiterm([nlmis  cR  1 -L(i)],    -co*l1'*C(:,:,l)',    c4*Om2);%%%%
    lmiterm([nlmis  cR  1 -L(i)],     co*l2'*C(:,:,l)',    c3*Om2);
    lmiterm([nlmis  cR  1 -L(i)],    -co*l2'*C(:,:,l)',    c4*Om2);
    lmiterm([nlmis  cR  1 -L(i)],     co*l3'*D(:,:,l)',    c3*Om2);
    lmiterm([nlmis  cR  1 -L(i)],    -co*l3'*D(:,:,l)',    c4*Om2);%%%%    
    
    lmiterm([nlmis  cR  1  W(i)],     co,                      Om2);
    %% (2,2)
    lmiterm([nlmis  cR  cC U(l,i)],   co,                      1);
    end
end
        
end