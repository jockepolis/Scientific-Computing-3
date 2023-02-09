% Our own LU-solver function

% Set the starting conditions
function myoLU = myownLU(A, b)
n = size(A,1);
% First iterate over the elements in the matrix A
for j = 1:n, i = j+1:n;
    A(i,j) = A(i,j)/A(j,j);
    A(i,i) = A(i,i)-A(i,j)*A(j,i);
end
% Now create the matrixes L and U by using matlab's tril and triu
L = tril(A,-1)+eye(n); % eye is the identity matrix since the -1 cancels out the diagonal
U = triu(A);

% Solve it as the formula; Ly = b and Ux = y
y = L\b; 
x = U\y;
myoLU = x;
end

