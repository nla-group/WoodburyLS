function x = QRUpdateLS(Q,R,U,V,b)
% QRUpdateLS Economy QR update to least squares solution .
% x = QRUpdateLS (Q,R,U,V,b) solves the least squares problem
% min_x || ( QR + UV ' )x - b ||, where Q is m x n orthogonal ,
% m >= n, R is n x n upper triangular , U is m x r, and V is n x r.
[n,r] = size(V); U0 = U;
% orthogonalise U wrt Q0
U = U-Q*(Q'*U);
[Qappend, ~] = qr(U ,0);
QTb = [Q'*b; Qappend'*b];
R = [R;zeros(r,n)];
W = [Q, Qappend]'*U0;
% reduce W to upper triangular using Givens rotations, and apply to R and QTb
for i = 1:r
    for k = n -1+i: -1:0+i
        [c,s] = Givens(W(k,i),W(k+1,i));
        G = [c s;-s c];
        idx = [k,k+1];
        R(idx ,1:n) = G'*R(idx ,1:n);
        W(idx ,:) = G'*W(idx ,:);
        QTb(idx) = G'* QTb(idx);
    end
end
R = R + W*V'; % R now has r subdiagonals
% now return R to upper triangular form using Householder reflectors
for i = 1:n
    [v,beta] = House(R(i:i+r,i));
    QTb(i:i+r) = QTb(i:i+r) - beta*v*(v'* QTb(i:i+r));
    R(i:i+r,i:n) = R(i:i+r,i:n) - (beta*v)*(v'*R(i:i+r,i:n));
end
R = R(1:n ,1:n);
x = R\QTb(1:n);
end

function [c,s] = Givens(a,b)
% function [c,s] = Givens(a,b)
% Givens rotation computation
% Determines cosine-sine pair (c,s) so that [c s;-s c]'*[a;b] = [r;0]
% GVL4: Algorithm 5.1.3
if b==0
    c = 1; s = 0;
else
    if abs(b)>abs(a)
        tau = -a/b; s = 1/sqrt(1+tau^2); c = s*tau;
    else
        tau = -b/a; c = 1/sqrt(1+tau^2); s = c*tau;
    end
end
end

function [v,beta] = House(x)
% function [v,beta] = house(x)
% Householder vector computation
% x is a column m-vector.
% v is a column m-vector with v(1) = 1 and beta is a scalar such
% that P = I - beta*v*v' is orthogonal and Px = norm(x)*e1.
% GVL4: Algorithm 5.1.1
m = length(x);
if m>1
    sigma = x(2:m)'*x(2:m); v = [1;x(2:m)];
    if sigma==0 && x(1)>0
        beta = 0;
    elseif sigma==0 && x(1)<0
        beta = 2;
    else
        mu = sqrt(sigma + x(1)^2);
        if x(1)<0
            v(1) = x(1) - mu;
        else
            v(1) = -sigma/(x(1)+mu);
        end
        beta = 2*v(1)^2/(sigma+v(1)^2);
        v = v/v(1);
    end
else
    v = 0; beta = 0;
end
end