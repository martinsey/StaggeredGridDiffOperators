clear all;

m = 5;
n = 5;
[grids, dx, dy, div, grad_div_, laplace] = DiffOperators(m, n);

% test calculation

p1 = [zeros(1, 5);randn(4, 5); zeros(1,5)];
p2 = [zeros(5, 1),randn(5, 4), zeros(5,1)];


dxp1 = p1(2:6, 1:5) - p1(1:5, 1:5);
dxp2 =  p2(1:5, 2:6) - p2(1:5, 1:5);
div = dxp1 + dxp2;
dxdiv = (div(2:5, 1:5) - div(1:4, 1:5)) * 25;
dydiv = (div(1:5, 2:5) - div(1:5, 1:4)) * 25;

result = grad_div_ * [p1(:); p2(:)];

assert(norm(reshape(result(1:30), 6, 5) - [zeros(1, 5);dxdiv;zeros(1, 5)]) < 1e-10);
assert(norm(reshape(result(31:60), 5, 6) - [zeros(5, 1),dydiv,zeros(5, 1)]) < 1e-10);

% test duality
m = 3; n = 3;

[grids, dx, dy, div_, grad_div_, laplace] = DiffOperators(m, n);
u = randn(3, 3);
s1 = randn(2, 3);
s2 = randn(3, 2);

s1_c = [zeros(1, 3);s1; zeros(1,3)];
s2_c = [zeros(3, 1),s2, zeros(3,1)];

dot([dx; dy] * u(:),[s1(:); s2(:)]);
dot(u(:), -div_ * [s1_c(:); s2_c(:)]);

assert(abs(dot([dx; dy] * u(:),[s1(:); s2(:)]) - dot(u(:), -div_ * [s1_c(:); s2_c(:)])) < 1e-10)
