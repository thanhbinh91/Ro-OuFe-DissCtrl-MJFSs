% clear
[nx, nu, nw, ny, nz, r, A, B, E, C, D1, G, H, D2] = SysParas;

T  = 100;
Ts = 2;
N  = T/Ts + 1;

Pi = 1;

x  = zeros(nx,N); x(:,1)  = [0.1;-0.1;0.1];
u  = zeros(nu,N);
xh = zeros(nx,N); xh(:,1) = [0; 0; 0];
y  = zeros(ny,N);
z  = zeros(nz,N);
th1  = zeros(1,N);
th1h = zeros(1,N);
th2  = zeros(1,N);
th2h = zeros(1,N);


mode = 1;
beta = 3.18;
[X, X_, P, F, L, copt] = SLPMM(0, beta);

% if mode == 1
%     Q = -eye(nz); Q1 = chol(-Q); nq = nz;
%     S =  zeros(nw,nz);
%     R = (beta+beta^2)*eye(nw);
% elseif mode == 2
%     Q = 0*eye(nz); Q1 = 0*eye(nz); nq = nz;
%     S = eye(nw,nz);
% else
%     Q = -0.36*eye(nz); Q1 = chol(-Q); nq = nz;
%     S = -2*ones(nw,nz);
%     R =  2*eye(nw);
% end

de = 0.01/pi;

Hinf = zeros(1,N);
Zsum = zeros(1,N);
Wsum = zeros(1,N);


for k = 1:N-1

    eta = x(2,k) - Ts/(2*5.5)*x(1,k);
    
%     if  k <= 10
%         wk = 1;
%     else
%         wk = 0;
%     end
    wk = exp(-0.3*k)*sin(2*k);
%     wk = 0; 
    %% Fuzzy basic functinos
    if eta~= 0
    th1(k) = (sin(eta) - de*eta)/(eta*(1-de));
    else
    th1(k) = 1;
    end
    th2(k) = 1 - th1(k);
      
    
    %% System matrices and output
    Ath  =  th1(k)*A(:,:,1)  +  th2(k)*A(:,:,2);
    Bth  =  th1(k)*B(:,:,1)  +  th2(k)*B(:,:,2);
    Eth  =  th1(k)*E(:,:,1)  +  th2(k)*E(:,:,2);
    Cth  =  th1(k)*C(:,:,1)  +  th2(k)*C(:,:,2);
    D1th =  th1(k)*D1(:,:,1) +  th2(k)*D1(:,:,2);
    D2th =  th1(k)*D2(:,:,1) +  th2(k)*D2(:,:,2);
    Gth  =  th1(k)*G(:,:,1)  +  th2(k)*G(:,:,2);
    Hth  =  th1(k)*H(:,:,1)  +  th2(k)*H(:,:,2);
    
    %% Mismatch fuzzy basic function
    if eta~= 0
    th1h(k) = (sin(eta) - de*eta)/(eta*(1-de));
    else
    th1h(k) = 1;
    end
    th2h(k) =  1 - th1h(k);
    
      
    %% Feedback and Observer Gains and Matrices
    Athh  =  th1h(k)*A(:,:,1)  +  th2h(k)*A(:,:,2);
    Bthh  =  th1h(k)*B(:,:,1)  +  th2h(k)*B(:,:,2);
    Ethh  =  th1h(k)*E(:,:,1)  +  th2h(k)*E(:,:,2);
    Cthh  =  th1h(k)*C(:,:,1)  +  th2h(k)*C(:,:,2);
    D1thh =  th1h(k)*D1(:,:,1) +  th2h(k)*D1(:,:,2);
    D2thh =  th1h(k)*D2(:,:,1) +  th2h(k)*D2(:,:,2);
    Gthh  =  th1h(k)*G(:,:,1)  +  th2h(k)*G(:,:,2);
    Hthh  =  th1h(k)*H(:,:,1)  +  th2h(k)*H(:,:,2);
    
    
    Pthh  =  th1h(k)*P(:,:,1)  +  th2h(k)*P(:,:,2); 
    Fthh  =  th1h(k)*F(:,:,1)  +  th2h(k)*F(:,:,2);
    Lthh  =  th1h(k)*L(:,:,1)  +  th2h(k)*L(:,:,2);
    
    
%     %%% Check LMIs solver %%%%
%     LMI21 = [-S*Gth - S*Hth*Fthh,  -S*Gth;
%              Q1*Gth + Q1*Hth*Fthh,  Q1*Gth];
%          
%     LMI31 = [Athh      +  Bthh*Fthh      + Lthh*(Cth-Cthh),   Lthh*Cth;
%             (Ath-Athh) + (Bth-Bthh)*Fthh - Lthh*(Cth-Cthh),   Ath-Lthh*Cth];
%         
%     LMI22 = [-S*D2th-D2th'*S'+beta*eye(nw)-R,  D2th'*Q1';
%               Q1*D2th,                        -eye(nz)];
%           
%     LMI32 = [Lthh*D1th,     zeros(nx,nz);
%              Eth-Lthh*D1th, zeros(nx,nz)];
%          
%     LMI  = [-Pthh,  LMI21', LMI31';
%              LMI21, LMI22,  LMI32';
%              LMI31, LMI32, -X(:,:,1)];
%          
%     if max(eig(LMI)) > 0
%         disp('Infeasible')
%     else
%     end
    
    u(k)   =  Fthh*xh(:,k);
    %% Observer
    y(:,k)   = Cth*x(:,k) + D1th*wk;
    xh(:,k+1) = Athh*xh(:,k) + Bthh*u(k) + Lthh*(y(:,k) - Cthh*xh(:,k));
    
    %% System
    x(:,k+1)  = Ath*x(:,k) + Bth*u(k) + Eth*wk;
    z(k)      = Gth*x(:,k) + Hth*u(k) + D2th*wk;
    
    Zsum(k+1) = Zsum(k) + z(k)'*z(k);
    Wsum(k+1) = Wsum(k) + wk'*wk;
    Hinf(k+1) = sqrt(Zsum(k+1)/(beta^2*Wsum(k+1)));
end

fontsize = 18;
linewidth = 1.5;

f1 = clf(figure(1)); axes('Position',[0.1 0.1 0.85 0.85]);
stairs(0:N-1,x(1,:),'linewidth',linewidth); hold on;
stairs(0:N-1,xh(1,:),'linewidth',linewidth, 'LineStyle','--');
set(gca,'fontsize',fontsize);
axis([0 N-1 -0.6 0.4]);
legend('$x_{1,k}$','$\hat x_{1,k}$','fontsize',fontsize,'interpreter','latex')
grid on;


clf(figure(2));  axes('Position',[0.1 0.1 0.85 0.85]);
stairs(0:N-1,x(2,:),'linewidth',linewidth); hold on;
stairs(0:N-1,xh(2,:),'linewidth',linewidth, 'LineStyle','--');
set(gca,'fontsize',fontsize);
axis([0 N-1 -0.2 0.22]);
legend('$x_{2,k}$','$\hat x_{2,k}$','fontsize',fontsize,'interpreter','latex')
grid on;


clf(figure(3)); axes('Position',[0.1 0.1 0.85 0.85]);
stairs(0:N-1,x(3,:),'linewidth',linewidth); hold on;
stairs(0:N-1,xh(3,:),'linewidth',linewidth, 'LineStyle','--');
set(gca,'fontsize',fontsize);
axis([0 N-1 -0.7 0.8]);
legend('$x_{3,k}$','$\hat x_{3,k}$','fontsize',fontsize,'interpreter','latex')
grid on;

clf(figure(4)); axes('Position',[0.1 0.1 0.85 0.85]);
stairs(0:N-1,u,'linewidth',linewidth);
legend('$u_k$','fontsize',fontsize,'interpreter','latex');
set(gca,'fontsize',fontsize);
axis([0 N-1 -0.8 1.1]);
grid on;


% clf(figure(5));
% stairs(0:N-1,th1-th1h,'linewidth',linewidth);
% set(gca,'fontsize',fontsize);
% legend('$|\theta_i(k) - \tilde \theta_i(k)|$','fontsize',fontsize,'interpreter','latex');

clf(figure(5)); axes('Position',[0.1 0.1 0.85 0.85]);
stairs(0:N-1,Hinf,'linewidth',linewidth);
set(gca,'fontsize',fontsize);
axis([0 N-1 0 1]);
legend('$\frac{\sum_{k=0}^{T}\Vert z_k \Vert_2^2}{\gamma_{min}^2 \sum_{k=0}^{T}\Vert w_k \Vert_2^2}$'...
    ,'fontsize',fontsize+2,'interpreter','latex','location','best');


