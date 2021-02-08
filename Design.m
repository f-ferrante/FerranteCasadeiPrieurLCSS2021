%This code solves the optimization problem considered in the paper. 
%The code requires YALMIP parser for Linear Matrix Inequality, freely avaialbe at https://yalmip.github.io. 
%Any SDP solver can be used.    

clear all;
%Plant data
L=[0 0 0 0; -1 2 -1 0; 0 -1 1 0; -1 0 0 1];%Laplacian of the network
%Convective speeds
lamb1=2;
lamb2=sqrt(2);
Lambda=blkdiag(lamb1,lamb2);
B=[2 0]';
H=[0 1;-1 0]; 
E=[.1;0];
Q=[1 0];
ne=min(size(E));
nq=min(size(Q));
mu=0.1474;
%smallest and largest positive eigenvalued of L
uLambda=0.3820;
oLambda=2.6180;

%%%%G=[0 S2; S2 -2I]
G12=1; %S2=1 for tanh
G22=-2*eye(ne);

%SDP variables
W=diag(sdpvar(2,1));
sigma=sdpvar(1,1);
Y=sdpvar(1,2);

%Conditions

Mat1=[-W*exp(-mu), (W*H'+uLambda*Y'*B');
    (W*H'+uLambda*Y'*B')',-W]; 

Mat2=[-W*exp(-mu), (W*H'+oLambda*Y'*B');
    (W*H'+oLambda*Y'*B')',-W];

Qs0=[-mu*W, sigma*inv(Lambda)*E+W*Q'*G12;
    (sigma*inv(Lambda)*E+W*Q'*G12)',sigma*G22
    ];

Qs1=[-mu*W*exp(-mu), sigma*exp(-mu)*inv(Lambda)*E+W*Q'*G12;
    (sigma*exp(-mu)*inv(Lambda)*E+W*Q'*G12)',sigma*G22
    ];

problem=[Mat1<=0*eye(max(size(Mat1))), Mat2<=0*eye(max(size(Mat2))), W>=1e-5*eye(2),Qs0<=-1e-8*eye(max(size(Qs0))),...
    Qs1<=-1e-8*eye(max(size(Qs1)))
    ];

%Solution to the LMI problem
options=sdpsettings('solver','sdpt3','verbose',3);
solution=solvesdp(problem,0,options);

K=double(Y)*inv(double(W)); %coupling gain
P=inv(double(W));           %matrix P defining the Lyapunov functional 
Thorizon=20;                %Time horizon for the simulation

global K P lamb1 lamb2 L Thorizon H B Q E 





