function kittekatt = ownLUtest(A,b)
%%%%%%%% Inputs %%%%%%%%
% A : state matrix
% B : input vector
%%%%%%%% Outputs %%%%%%%%
% x_soln : solution vector to U*x_soln=xstar_soln
% xstar_soln: solution vector to L*xstar_soln=b
% L : Lower triangular matrix
% U : Upper triangular matrix
n=rank(A);
%Initialize L and U matrix
L=zeros(n);
U=eye(n);
for s=1:n
%Perform calculation on column j to determine L components 
j=s;
for i=j:n
    L(i,j)=A(i,j)-L(i,1:(j-1))*U(1:(j-1),j);
end
%Perform calculation on row i to determine U components   
i=s;
U(i,i)=1;
for j=i+1:n
    U(i,j)=(A(i,j)-L(i,1:(i-1))*U(1:(i-1),j))/(L(i,i));
end
end
%obtain xstar_soln by solving L*xstar_soln=A using forward substitution
xstar_soln(n)=0;
xstar_soln=xstar_soln';
xstar_soln(1)=b(1)/L(1,1);
for i=2:n
    xstar_soln(i)=(b(i,1)-L(i,1:i-1)*xstar_soln(1:i-1))/L(i,i);
end
%obtain x_soln by solving U*x_soln=xstar_soln using backward substitution
x_soln(n)=0;
x_soln=x_soln';
x_soln(n)=xstar_soln(n);
i=n-1;
while i>0
    x_soln(i)=xstar_soln(i)-U(i,i+1:n)*x_soln(i+1:n);
   i=i-1;
end
kittekatt = x_soln;