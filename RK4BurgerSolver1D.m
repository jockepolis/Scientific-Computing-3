function  RK4BurgerSolver1D()
    clear all; close all; clc;
    n=input(' Choose an N = {41, 81, 161, 321, 641} : ');
    a = -1; b = 1; c = 2; eps = 0.1;
    h = abs(a-b)/(n-1); % Spatial steps
    m = 5000; % Number of time levels
    T = 0.4; %Final time
    t = linspace(0,T,m+1); % timevector
    x = linspace(a,b,n); % Space vector
    f = @(x,t,eps) (c-tanh((x + 1/2 - c.*t)./(2*eps))); %Exact function
    BETA = @(u) (1/2.* u .^2);
    u_0 = f(x,0,eps); %Initial data
    
    M = Mass_assembler(x);
    A = Advection_assembler(x);
    S = StiffnessAssembler1D(x);
    
    U=zeros(length(x),length(t));
    U(:,1) = u_0;
    k = t(2) - t(1);
    for i = 1:m
        K1 = cg(M,(-A*BETA(U(:,i)) - eps*S*U(:,i)), 1e-3);
        
        K2 = cg(M,(-A*BETA(U(:,i)+k/2*K1) - eps*S*(U(:,i)+k/2*K1)), 1e-3);
        
        K3 = cg(M,(-A*BETA(U(:,i)+k/2*K2) - eps*S*(U(:,i)+k/2*K2)),1e-3);
        
        K4 = cg(M,(-A* BETA((U(:,i)+k*K3)) - eps*S*(U(:,i)+k*K3)), 1e-3);
        U(:,i+1)=U(:,i) + k/6*(K1+2*K2+2*K3+K4);
        
        U(1,i+1) = f(a,t(i+1),eps);
        U(end,i+1) = f(b,t(i+1),eps);
    end
    mesh(x,t,U');
    xlabel('Space'); ylabel('Time'), title('Solution')
end