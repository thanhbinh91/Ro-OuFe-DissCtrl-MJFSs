% clear
% clear
[nx, nu, nw, ny, nz, s, r, A, B, E, C, D, G, H, Jj] = SysParas;

T  = 30;
Ts = 0.5;
N  = T/Ts + 1;
Ns = 1:1:s;

Pi = [0.8 0.1 0.1;
      0.2 0.7 0.1;
      0.5 0.2 0.3];

x  = zeros(nx,N); x(:,1) = [0.2*pi;-0.5];
u  = zeros(nu,N);
xh = zeros(nx,N); xh(:,1) = [0; 0];
y  = zeros(ny,N);
z  = zeros(nz,N);
mode = ones(1,N);
th1  = zeros(1,N);
th1h = zeros(1,N);
th2  = zeros(1,N);
th2h = zeros(1,N);

beta = 3.65;
%[X, X_, P, F, L, topt] = SLPMM(0.2, beta);

de = 0.01/pi;


for k = 1:N-1
    g = mode(k);
%     if k >= 55 && k <= 60
%     wk1 = 5;
%     else
%     wk1 = 0;
%     end
%     wk2 = exp(-3*k*Ts)*sin(5*k*Ts);
%     wk  = [wk1;wk2];
    wk = 0.4*exp(-0.1*k)*sin(0.3*k);
     
    %% Fuzzy basic functinos
    if x(1,k) ~= 0
        th1(k) = (sin(x(1,k)) - de*x(1,k))/((1-de)*x(1,k));
    else
        th1(k) = 1;
    end
    th2(k) = 1 - th1(k);
      
    
    %% System matrices and output
    Agth  =  th1(k)*A(:,:,g,1)  +  th2(k)*A(:,:,g,2);
    Bgth  =  th1(k)*B(:,:,g,1)  +  th2(k)*B(:,:,g,2);
    Egth  =  th1(k)*E(:,:,g,1)  +  th2(k)*E(:,:,g,2);
    Cgth  =  th1(k)*C(:,:,g,1)  +  th2(k)*C(:,:,g,2);
    Dgth  =  th1(k)*D(:,:,g,1)  +  th2(k)*D(:,:,g,2);
    Ggth  =  th1(k)*G(:,:,g,1)  +  th2(k)*G(:,:,g,2);
    Hgth  =  th1(k)*H(:,:,g,1)  +  th2(k)*H(:,:,g,2);
    
    %% Mismatch fuzzy basic function
    y(k) = Cgth*x(:,k) + Dgth*wk;
    if xh(1,k) ~= 0
        th1h(k) = (sin(xh(1,k)) - de*xh(1,k))/((1-de)*xh(1,k));
    else
        th1h(k) = 1;
    end
    th2h(k) = 1 - th1h(k);  
    
      
    %% Feedback and Observer Gains and Matrices
    Agthh =  th1h(k)*A(:,:,g,1)  +  th2h(k)*A(:,:,g,2);
    Bgthh =  th1h(k)*B(:,:,g,1)  +  th2h(k)*B(:,:,g,2);
    Egthh =  th1h(k)*E(:,:,g,1)  +  th2h(k)*E(:,:,g,2);
    Cgthh =  th1h(k)*C(:,:,g,1)  +  th2h(k)*C(:,:,g,2);
    Dgthh =  th1h(k)*D(:,:,g,1)  +  th2h(k)*D(:,:,g,2);
    Ggthh =  th1h(k)*G(:,:,g,1)  +  th2h(k)*G(:,:,g,2);
    Hgthh =  th1h(k)*H(:,:,g,1)  +  th2h(k)*H(:,:,g,2);
     
    Fgthh =  th1h(k)*F(:,:,g,1)  +  th2h(k)*F(:,:,g,2);
    Lgthh =  th1h(k)*L(:,:,g,1)  +  th2h(k)*L(:,:,g,2);
    
    u(k) = Fgthh*xh(:,k);
    %% Observer
    y(k) = Cgth*x(:,k) + Dgth*wk;
    xh(:,k+1) = Agthh*xh(:,k) + Bgthh*u(k) + Lgthh*(y(k) - Cgthh*xh(:,k));
    
    %% System
    x(:,k+1)  = Agth*x(:,k) + Bgthh*u(k) + Egthh*wk;
    
    %% Markov
    if mod(k,1) == 0
        mode(k+1) = randsample(Ns,1,true,Pi(g,:));
    else
        mode(k+1) = mode(k);
    end
end

fontsize = 16;
linewidth = 1.5;

% clf(figure(16)); axes('Position',[0.1 0.1 0.85 0.85]);
% f1 = stairs(0:N-1,mode,'linewidth',linewidth,'color','k');
% yticks([1 2 3]); axis([0 N-1 0.8 3.4]);
% set(gca,'fontsize',fontsize);
% legend('$\phi_k$','fontsize',fontsize,'interpreter','latex')


clf(figure(17)); axes('Position',[0.1 0.1 0.85 0.85]);
stairs(0:N-1,x(1,:),'linewidth',linewidth); hold on;
stairs(0:N-1,xh(1,:),'linewidth',linewidth,'LineStyle','--');
set(gca,'fontsize',fontsize);
axis([0 N-1 -0.2 0.9]);
legend('$x_{1,k}$','$\hat x_{1,k}$','fontsize',fontsize,'interpreter','latex');
grid on


clf(figure(18)); axes('Position',[0.1 0.1 0.85 0.85]);
stairs(0:N-1,x(2,:),'linewidth',linewidth); hold on;
stairs(0:N-1,xh(2,:),'linewidth',linewidth,'LineStyle','--');
set(gca,'fontsize',fontsize);
axis([0 N-1 -1.1 0.5]);
legend('$x_{2,k}$','$\hat x_{2,k}$','fontsize',fontsize,'interpreter','latex')
grid on


clf(figure(19)); axes('Position',[0.1 0.12 0.85 0.85]);
stairs(0:N-1,u,'linewidth',linewidth);
legend('$u_k$','fontsize',fontsize,'interpreter','latex');
set(gca,'fontsize',fontsize);
grid on


clf(figure(20)); axes('Position',[0.1 0.1 0.85 0.85]);
stairs(0:N-1,abs(th1-th1h),'linewidth',linewidth);
set(gca,'fontsize',fontsize);
axis([0 N-1 -0.01 0.05]);
legend('$|\theta_{1,k} - \hat{\theta}_{1,k}|$','fontsize',fontsize+2,'interpreter','latex','location','best');
grid on


clf(figure(21)); axes('Position',[0.1 0.12 0.85 0.85]);
stairs(0:N-1,mode,'linewidth',linewidth);
legend('$\rho_k$','fontsize',fontsize,'interpreter','latex');
set(gca,'fontsize',fontsize);
set(gca,'YTick',[1 2 3],'YLim',[0.8 3.3]);
grid on



