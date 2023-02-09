clear all; close all; clc;
g=0;
L2NORM=zeros(5,1);
T=1;
for N=[41, 81, 161, 321, 641]
    g=g+1;
    x=linspace(0,1,N);
    h=1/(N-1);
    k=0.1*h;
    u=zeros(N,round(T/k));
    A=diag(-ones(N-1,1),-1)+diag(ones(N-1,1),1);
    A(1,end)=-1;
    A(end,1)=1;
    A=1/(2*h).*A;
    u_0=zeros(N,1);
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
    
    u(:,1)=u_0;
    for i=1:T/k
        K1=PDE(u(:,i),A);
        K2=PDE(u(:,i)+k/2*K1,A);
        K3=PDE(u(:,i)+k/2*K2,A);
        K4=PDE(u(:,i)+k*K3,A);
        u(:,i+1)=u(:,i)+k/6*(K1+2*K2+2*K3+K4);
    end
L2NORM(g)=norm(((u(:,end)'-u(:,1)')./sqrt(N)));
h_plot(g)=h;
end
%%
loglog(1./([41 81 161 321 641]-1),L2NORM,'ok-','MarkerFaceColor',[0.635 0.078 0.1840]);
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex')
xlabel('Spatial step h','Interpreter','latex');ylabel('$|| u-u_h ||_2$','Interpreter','latex');
title('$ || u-u_h ||_2$ plotted against the spatial step h','Interpreter','latex');
grid on
axis padded
%%
q=log(L2NORM(end)/L2NORM(end-1))/log(h_plot(end)/h_plot(end-1))

%% Functions
function dudt=PDE(u,A)
dudt=-1.*A*u;
end