function [omega1, omega2] = grid(m ,n)

% implementation for non-square grids follow later
assert(m == n);

[X1, Y1] = ndgrid(0:m, 1:n);
[X2, Y2] = ndgrid(1:m, 0:n);

XX1 = X1(:);
YY1 = Y1(:);

XX2 = X2(:);
YY2 = Y2(:);

omega1 = {XX1, YY1};
omega2 = {XX2, YY2};