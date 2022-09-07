function [X, X_, P, F, L, copt] = SLPMM(al, beta)
mode = 3;
Nmax = 100;

[nx, nu, nw, ny, nz, s, r, A, B, E, C, D, G, H, J, Pi] = SysParas;


%% Step 1
[Xk, Xk_] = Theorem3_S1(mode, al, beta);
disp('=== Theorem 3: Step 2 LMIs solvers ===')
L  = 0;
F  = 0;
P  = 0;
X  = 0;
X_ = 0;
copt = 0;
for k = 1:Nmax
    [Xo, Xo_, Po, Fo, Lo, to] = Theorem3_S2(mode, al, beta, Xk, Xk_);
    fprintf('+ minTrace: t%d = %4.3f,',k,to);

    copt = to;  

    if copt - 4*nx < 0.01
        L    = Lo;
        F    = Fo;
        P    = Po;
        X    = Xo;
        X_   = Xo_;
        
        if mode == 1
            fprintf('\n Feasible in H_inf performance with %4.3f\n', beta);
        elseif mode == 2
            fprintf('\n Feasible in Passivity performance with %4.3f\n', beta);
        else
            Q = -0.01*eye(nz); Q1 = chol(-Q); nq = nz;
            S =  0.2*ones(nw,nz);
            R =  5*eye(nw);
            fprintf(['\n FEASIBLE in (Q=%4.2f,S=%4.2f,R=%4.2f)-' ...
                'dissipative performance with %4.3f\n'], Q, S, R, beta);
        end
        
        break
    end
    
    a = 0;
    b = 0;
    c = 0;
    
    a = a + trace((Xo - Xk)*(Xo_- Xk_));
    b = b + trace((Xo - Xk)*Xk_ + Xk*(Xo_ - Xk_));
    c = c + trace(Xk*Xk_);
    
    f  = @(x) a*x^2+b*x+c;
    ga = fminbnd(f, 0, 1);
    fprintf('   gamma = %4.3f\n',ga);
    if ga < 0.001
        disp('The set of LMIs is INFEASIBLE');
        break
    end
    Xk  = (1-ga)*Xk  + ga*Xo;
    Xk_ = (1-ga)*Xk_ + ga*Xo_;
end
fprintf('\n');
end