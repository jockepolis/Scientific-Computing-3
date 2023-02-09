%% INITIAL VALUES ETC.
clear all; close all; clc;

n=input(' Choose an N = {41, 81, 161, 321, 641} : ');
a = -1;
b = 1;
c = 2;
epsilon = 0.1;
h=abs(a-b)/(n-1);
T=0.4;
k=0.1*h;
t=0:k:T;
x=linspace(a,b,n);
A=AdvectionAssembler(x);
S=my_stiffness_matrix_assembler(x);
M=MassAssembler1D(x);

u_0=c-tanh((x+1/2)/(2*epsilon));

u_exact=zeros(n,round(T/k));
i=0;
for tt=linspace(0,T,round(T/k)+1)
    i=i+1;
    u_exact(:,i)=c-tanh((x+1/2-c.*tt)/(2*epsilon));
end

g_a=u_exact(1,:); g_b=u_exact(end,:);
g_a_x=-(sech(a+1/2-c.*t)/(2*epsilon)).^2/(2*epsilon); %Derivative of g_a with respect to x
g_b_x=-(sech(b+1/2-c.*t)/(2*epsilon)).^2/(2*epsilon);

hat_0=zeros(n,1)'; hat_n=zeros(n,1)';
hat_0(1)=1; hat_n(end)=1;
b_vec=zeros(n,1); c_vec=zeros(n,1); d_vec=zeros(n,1);

%% RK4
u=zeros(n,round(T/k));
u(:,1)=u_0;
for i=1:round(T/k)
    b_vec(1)=g_a_x(i); b_vec(end)=g_b_x(i);
    c_vec(1)=BETA(g_a(i)); c_vec(end)=BETA(g_b(i));
    d_vec(1)=g_a(i) ; d_vec(end)=g_b(i);
    
    K1=PDE(u(:,i),A,@BETA,epsilon,S,b_vec,c_vec,d_vec); %PDE(u,A,BETA,epsilon,S,b,c,d)
    K2=PDE(u(:,i)+k/2*K1,A,@BETA,epsilon,S,b_vec,c_vec,d_vec);
    K3=PDE(u(:,i)+k/2*K2,A,@BETA,epsilon,S,b_vec,c_vec,d_vec);
    K4=PDE(u(:,i)+k*K3,A,@BETA,epsilon,S,b_vec,c_vec,d_vec);
    bb = k/6*(K1+2*K2+2*K3+K4);
    u(:,i+1)=u(:,i)+M\bb;
    
%     plot(x,u(:,i+1));
%     title(['Solution at time= ',num2str(round(t(i+1),2))]);
%     drawnow limitrate
%     F(i+1)=getframe(gcf);
end

%% PLOTTING
mesh(x,t,u');
xlabel('Space'); ylabel('Time'), title('Solution')

%% FUNCTIONS
function dudt=PDE(u,A,BETA,epsilon,S,b,c,d)
dudt=(-A*BETA(u)-epsilon.*S*u-b-c-epsilon.*d);
end

function M = MassAssembler1D(x)
n = length(x)-1;
M = zeros(n+1,n+1);
for i = 1:n 
    h = x(i+1) -x(i);
    M(i,i) = M(i,i) + h/3;
    M(i,i+1) = M(i,i+1) + h/6;
    M(i+1,i) = M(i+1,i) + h/6;
    M(i+1,i+1) = M(i+1,i+1) + h/3;
end
end

% function d = LoadAssembler1D(x,f)
% n = length(x)-1;
% d = zeros(n+1,1);
% for i = 1:n
%     h = x(i+1) -x(i);
%     d(i) = d(i) + f(x(i))*h/2;
%     d(i+1) = d(i+1) + f(x(i+1))*h/2;
% end
% end

function A = AdvectionAssembler(x)
n = length(x)-1;
A = (zeros(n+1,n+1)+diag(-0.5*ones(n, 1),-1)+diag(0.5*ones(n,1),1));
end

function S=my_stiffness_matrix_assembler(x)
n = length(x) - 1; % number of elements
S = zeros(n+1, n+1); % initialize stiffnes matrix to zero
for i = 1:n % loop over elements
    h = x(i+1) - x(i); % element length
    n = [i i+1]; % nodes
    S(n,n) = S(n,n) + [1 -1; -1 1]/h; % assemble element stiffness
end
 S(1,1) = 1.e+6; % adjust for BC
 S(n,n) = 1.e+6;
end

function b = BETA(u)
b = 1/2*u.*u;
end
