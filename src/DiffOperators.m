function [grids, dx, dy, div, grad_div, laplace_d] = DiffOperators(m, n)
%LAPLACE Get staggered grid differential operators currently only
%implemented for square, i.e. m == n. 

[omega1, omega2] = Grid(m, n);
grids = {omega1, omega2};

% grad div
div = [sparse(m, 1),[speye(m - 1); sparse(1, m - 1)] - [sparse(1, m - 1); speye(m - 1)], sparse(m, 1)];
div = repmat({div}, 1, n);
div = [blkdiag(div{:}),sparse(m*n, (m + 1) * n)]  ...
    + [sparse(m*n, (m + 1) * n),[sparse(m * n, m), [speye((m - 1) * n); sparse(n, (m - 1) * n)], sparse(m * n, m)]...
    + [sparse(m * n, m), [sparse(n, (m - 1) * n); -speye((m - 1) * n)], sparse(m * n, m)]];

dx = repmat({-speye(m - 1, n) + [sparse(m - 1, 1), speye(n-1, m-1)]}, 1, m);
dx = blkdiag(dx{:});
dy = -speye((m - 1) * n, m * n) + [sparse((m - 1) * n, m), speye((m - 1) * n)];

grad_div = [dx; dy] * div;

div = div * m;

dx = dx * m;
dy = dy * n;

% insert zero rows for boundary nodes
bd = repmat({[sparse(1, m - 1);speye(m - 1);sparse(1, m - 1)]}, 1, m);
bd = blkdiag(bd{:});
grad_div1 = bd * grad_div(1:(m-1)*n,1:2*(m+1)*n);
grad_div2 = [sparse(n, 2*(m + 1)*n);grad_div((m - 1)*n + 1:2*(m - 1)*n,:); sparse(n, 2*(m + 1)*n)];

grad_div = [grad_div1;grad_div2] * (m * n); 

%laplace operator
laplace_block = repmat({[sparse(m + 1, 1),[sparse(1, m - 1);spdiags([ones(m - 1, 1) -4*ones(m - 1, 1) ones(m - 1, 1)], -1:1, m - 1, m - 1); sparse(1, m - 1)], sparse(m + 1, 1)]}, 1, m);
laplace_diag = blkdiag(laplace_block{:});
subdiag = repmat({[sparse(n + 1, 1), [sparse(1, n - 1); speye(m - 1); sparse(1, n - 1)], sparse(n + 1, 1)]},1, m - 1);
subdiag1 = [[sparse((m - 1)*(n + 1), m + 1), blkdiag(subdiag{:})]; sparse(n + 1, (m + 1)*n)];
subdiag2 = [sparse(n + 1, (m + 1)*n); [blkdiag(subdiag{:}), sparse((m - 1)*(n + 1), m + 1)]];
laplace1 = laplace_diag + subdiag1 + subdiag2;

diag = repmat({spdiags([ones(m, 1) -4*ones(m, 1) ones(m, 1)], -1:1, m, m)}, 1, m - 1);
laplace_diag = blkdiag(diag{:});

laplace_diag = spdiags([ones(n * (m - 1), 1), ones(n * (m - 1), 1)],[-n, n], n * (m - 1),m*(n - 1)) + laplace_diag;
laplace2 = [sparse((m + 1)*n, m), [sparse(m, (m - 1)*n); laplace_diag; sparse(m, (m - 1)*n)], sparse((m + 1)*n, m)];

laplace_d = [laplace1 * (m * n), sparse((m + 1) * n,  (n + 1) * m); sparse(m * (n + 1), (m + 1) * n), laplace2 * (m * n)];

end
