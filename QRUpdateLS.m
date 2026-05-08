function x = QRUpdateLS(Q,R,U,V,b)
%QRUpdateLS Economy QR update to least squares solution.
% x = QRUpdateLS (Q,R,U,V,b) solves the least squares problem
% min_x || ( QR + UV ' )x - b ||, where Q is m x n orthogonal,
% m >= n, R is n x n upper triangular , U is m x r, and V is n x r.
% We assume m >= n+r.
% Bloor Riley et al, A Sherman–Morrison–Woodbury Approach to Solving
% Least Squares Problems with Low-Rank Updates, 2026.

[n,r] = size(V);
% orthogonalise U wrt Q
U1 = U-Q*(Q'*U);
Qappend = orth(U1,10*eps);
QTb = [Q'*b; Qappend'*b];
R = [R; zeros(size(Qappend,2),n)];
W = [Q'*U; Qappend'*U];
%d = min(size(W,1), n+r);
% reduce W to upper triangular using Givens rotations, and apply to R and QTb
for i = 1:r
    for k = size(W,1)-1:-1:i
        idx = [k,k+1];
        G = planerot(W(idx,i));
        R(idx,:) = G*R(idx,:);
        W(idx,:) = G*W(idx,:);
        QTb(idx) = G*QTb(idx);
    end
end
R = R + W*V'; % R is now upper trapezoidal
% now return R to upper triangular form using Householder reflectors
for k = 1:(size(R,1)-r)
    [v,beta] = gallery("house", R(k:k+r,k));
    R(k:k+r,k:n) = R(k:k+r,k:n) - (beta*v)*(v'*R(k:k+r,k:n));
    QTb(k:k+r) = QTb(k:k+r) - beta*v*(v'*QTb(k:k+r));
end
R = R(1:n,1:n);
x = R\QTb(1:n);
end