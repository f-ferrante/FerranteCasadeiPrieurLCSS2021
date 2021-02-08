
function [x1,t1, V1to, V2to,V3to, V4to, V5to, V6to, V7to, V8to]=hyp_static_bc
global  K B H lamb1 lamb2 Thorizon
method = 'LxF';
%method = 'SLxW';
pdefun = @pdes;
Nx=500;
x = linspace(0,1,Nx);
Tinit = 0;
t = Tinit;

% Initial conditions %%%% SHOULD SATISFY THE COMPATIBILITY CONDITION
V(1,:)=(cos(2*pi*x)-1)*1;
V(2,:)=-(cos(2*pi*x)-1)*1;

V(3,:)=0*x.*(x-1);
V(4,:)=-0*(sin(2*pi*x));

V(5,:)=2*(sin(2*pi*x));
V(6,:)=-0*(sin(4*pi*x));

V(7,:)=0*(sin(2*pi*x));
V(8,:)=-2*(sin(4*pi*x));

sol = setup(1,pdefun,t,x,V,method,[],@bcfun);

% check CFL condition
dx = x(2)-x(1);
%dt = 0.95*dx/lamb2;
dt = 0.95*dx/2;
dt = 0.7*dx/2;
%dt = 0.3*dx/2;

Nt=floor(Thorizon/dt);
howfar = Thorizon/Nt;

% to stack the solution
V1to =zeros(Nt,Nx);
V2to =zeros(Nt,Nx);
V3to =zeros(Nt,Nx);
V4to =zeros(Nt,Nx);
V5to =zeros(Nt,Nx);
V6to =zeros(Nt,Nx);
V7to =zeros(Nt,Nx);
V8to =zeros(Nt,Nx);


V1to(1,:) = V(1,:);
V2to(1,:) = V(2,:);


V3to(1,:) = V(3,:);
V4to(1,:) = V(4,:);

V5to(1,:) = V(5,:);
V6to(1,:) = V(6,:);

V7to(1,:) = V(7,:);
V8to(1,:) = V(8,:);

for m = 1:Nt-1
    disp(['time step: ' num2str(m) '/' num2str(Nt)]) 
    %keyboard
    sol = hpde(sol,howfar,dt);  
    t = sol.t;
    V = sol.u;
    if 0
        figure(Nfig)
        plot(x,V,'.-k')
        axis(scale_axis) 
        title(['Solution at t = ',num2str(t),'.']);
        pause(1)
    end
    V1to(m+1,:) = V(1,:);
    V2to(m+1,:) = V(2,:);
    V3to(m+1,:) = V(3,:);
    V4to(m+1,:) = V(4,:);
    V5to(m+1,:) = V(5,:);
    V6to(m+1,:) = V(6,:);
    V7to(m+1,:) = V(7,:);
    V8to(m+1,:) = V(8,:);
end
t = linspace(Tinit,Thorizon,Nt);
[x1,t1] = meshgrid(x,t);
    
%=========================================================================
% Subfunctions

    function F = pdes(t,x,V,V_x) % linear PDE
         %global lamb1 lamb2
         %keyboard
        F = zeros(size(V));
       % f1=abs(V(1,:)+1)-abs(V(1,:)-1);
        %f2=abs(V(3,:)+1)-abs(V(3,:)-1);
        %f3=abs(V(5,:)+1)-abs(V(5,:)-1);
        f1=0.1*sin(V(1,:));
        f2=0.1*sin(V(3,:));
        f3=0.1*sin(V(5,:));
        f4=0.1*sin(V(7,:));
        
        %f1=0*tanh(V(1,:));
        %f2=0*tanh(V(3,:));
        %f3=0*tanh(V(5,:));
        f1=0.1*tanh(V(1,:));
        f2=0.1*tanh(V(3,:));
        f3=0.1*tanh(V(5,:));
        f4=0.1*tanh(V(7,:));
        
        F(1,:) = -lamb1.*V_x(1,:)+f1;
        F(2,:) = -lamb2.*V_x(2,:);
        
        F(3,:) = -lamb1.*V_x(3,:)+f2;
        F(4,:) = -lamb2.*V_x(4,:);
        
        F(5,:) = -lamb1.*V_x(5,:)+f3;
        F(6,:) = -lamb2.*V_x(6,:);
        
        F(7,:) = -lamb1.*V_x(7,:)+f4;
        F(8,:) = -lamb2.*V_x(8,:);
    % end function pdes
    end
    function [XL,XR] = bcfun(t,XLex,XRex)
    global  L
    XLtmp = XLex;
    %keyboard

    Xbc1=H*[XRex(1);XRex(2)]+B*K*(L(1,1)*[XRex(1);XRex(2)]...
         +L(1,2)*[XRex(3);XRex(4)]+L(1,3)*[XRex(5);XRex(6)]+L(1,4)*[XRex(7);XRex(8)]);
     
   
     Xbc2=H*[XRex(3);XRex(4)]+B*K*(L(2,1)*[XRex(1);XRex(2)]...
         +L(2,2)*[XRex(3);XRex(4)]+L(2,3)*[XRex(5);XRex(6)]+...
         L(2,4)*[XRex(7);XRex(8)]);
    
     Xbc3=H*[XRex(5);XRex(6)]+B*K*(L(3,1)*[XRex(1);XRex(2)]...
         +L(3,2)*[XRex(3);XRex(4)]+L(3,3)*[XRex(5);XRex(6)]+...
         L(3,4)*[XRex(7);XRex(8)]);
     
     Xbc4=H*[XRex(7);XRex(8)]+B*K*(L(4,1)*[XRex(1);XRex(2)]...
         +L(4,2)*[XRex(3);XRex(4)]+L(4,3)*[XRex(5);XRex(6)]+...
         L(4,4)*[XRex(7);XRex(8)]);
    
    XLtmp(1) = [1 0]*Xbc1;
    XLtmp(2) = [0 1]*Xbc1;
    XLtmp(3) = [1 0]*Xbc2;
    XLtmp(4) = [0 1]*Xbc2;
    XLtmp(5) = [1 0]*Xbc3;
    XLtmp(6) = [0 1]*Xbc3;
    XLtmp(7) = [1 0]*Xbc4;
    XLtmp(8) = [0 1]*Xbc4;
    XRtmp = XRex;
    XL=XLtmp;
    XR=XRtmp;
    end
end