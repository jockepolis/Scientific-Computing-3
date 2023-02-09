clear all; close all; clc;
g=0;
L2NORM=zeros(5,1);
T=1;
for N=[41, 81, 161, 321, 641, 1281]
    g=g+1;
    h=1/(N-1);
    k=0.5*h;
    t=0:k:T;
    x=linspace(0,1,N);
    u_0=zeros(1,N);
% With pen and paper, Epsilon was determined to be Epsilon = h/2
% for the discretization to become a upwind method
    epsilon = h/2;
    for i=1:N
        if abs(2*x(i)-0.3) <= 0.25
            u_0(i) = exp(-300*(2*x(i)-0.3)^2);
        elseif abs(2*x(i)-0.9) <= 0.2
            u_0(i) = 1;
        elseif abs(2*x(i)-1.6) <=0.2
            u_0(i) = sqrt(1-((2*x(i)-1.6)/0.2)^2);
        else
            u_0(i) = 0;
        end
    end

%% Central difference discretization in space
    u=zeros(N,100);
    u(:,1)=u_0;
    global A
    A=diag((1/(2*h)+epsilon/h^2).*ones(N-1,1),-1)+diag((-2*epsilon/(h^2)).*ones(N,1));
    A(1,end)=1/(2*h)+epsilon/h^2;
%% Forward Euler - discretization in time

    for j=1:T/k
        u(:,j+1)=u(:,j) + k*A*u(:,j);
    end
L2NORM(g)=norm(((u(:,end)'-u(:,1)')./sqrt(N)));
h_plot(g)=h;
end
%%
loglog(1./([41 81 161 321 641 1281]-1),L2NORM,'ok-','MarkerFaceColor',[0.635 0.078 0.1840]);
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex')
xlabel('Spatial step h','Interpreter','latex');ylabel('$|| u-u_h ||_2$','Interpreter','latex');
title('$ || u-u_h ||_2$ plotted against the spatial step h','Interpreter','latex');
grid on
axis padded
%%
q=log10(L2NORM(end)/L2NORM(end-1))/log10(h_plot(end)/h_plot(end-1))

%% Functions
