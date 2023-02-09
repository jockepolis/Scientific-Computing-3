function my_first_fem_solver()
a = -1; % left end point of interval
b = 1; % right
c = 2;
epsilon = 0.1;
N = 56; % number of intervals
h = (b-a)/N; % mesh size
x = a:h:b; % node coords
A=my_stiffness_matrix_assembler(x);
B=my_load_vector_assembler(x);
xi = A\B; % solve system of equations
plot(x,xi, 'r') % plot solution
hold on
plot(x,2-tanh((x+0.5-2*0.4)/(2*0.1)), 'b')
