function x = QRUpdateLS(Q,R,U,V,b)
% QRUpdateLS Economy QR update to least squares solution .
% x = QRUpdateLS (Q,R,U,V,b) solves the least squares problem
% min_x || ( QR + UV ' )x - b ||, where Q is m x n orthogonal ,
% m >= n, R is n x n upper triangular , U is m x r, and V is n x r.
[n,r] = size(V);
% orthogonalise U wrt Q0
U1 = U-Q*(Q'*U); 
U1 = U1-Q*(Q'*U1); % perform orthogonalisation twice
[Qappend, ~] = qr(U1 ,0);
QTb = [Q'*b; Qappend'*b];
R = [R;zeros(r,n)];
W = [Q, Qappend]'*U;
% reduce W to upper triangular using Givens rotations, and apply to R and QTb
for i = 1:r
    for k = n-1+r:-1:i
        idx = [k,k+1];
        G = planerot(W(idx,i));
        R(idx,1:n) = G*R(idx,1:n);
        W(idx,:) = G*W(idx,:);
        QTb(idx) = G*QTb(idx);
    end
end
R = R + W*V'; % R(1:n+r,1:n) is now upper trapezoidal
% now return R to upper triangular form using Householder reflectors
for i = 1:n
    [v,beta] = gallery("house", R(i:i+r,i));
    QTb(i:i+r) = QTb(i:i+r) - beta*v*(v'*QTb(i:i+r));
    R(i:i+r,i:n) = R(i:i+r,i:n) - (beta*v)*(v'*R(i:i+r,i:n));
end
R = R(1:n ,1:n);
x = R\QTb(1:n);
end