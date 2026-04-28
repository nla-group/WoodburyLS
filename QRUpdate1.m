function [Q,R] = QRUpdate1(Q0,R0,U,V)
%Rank-r update of a QR factorisation by applying the rank-1 update r times
% Q0 is mxn orthogonal, R0 is nxn upper triangular, U is mxr,
% and V is nxr.
% Q is mxn orthogonal and  R is nxn upper triangular so that
% QR = Q0*R0 + U*V'.
n = length(R0); r= size(U,2);
Q = Q0; R = R0;
for i = 1:r   %Iterate over r columns of U and V
    qperp = U(:,i) - Q(:,1:n)*(Q(:,1:n)'*U(:,i)); %Orthogonalise ith column of U
    Q(:,n+1) = qperp/norm(qperp);                 %normalises and appends to  Q
    R(n+1,:) = zeros(1,n);                        %Adds a row of zeros to R
    [Q,R] = Rank1Update(Q, R, U(:,i), V(:,i));    %Applies rank-1 update to QR
end
Q = Q(:,1:n); R = R(1:n,:);  %Removes added columns and rows of Q and R.
end


function [Q,R] = Rank1Update(Q0,R0,u,v)
% function [Q,R] = Rank1Update(Q0,R0,u,v)
% Rank-1 update of a QR factorization.
% Q0 is mxn orthogonal, R0 is mxn upper triangular, u is mx1,
% and v is nx1.
% Q is mxm orthogonal and  R is mxn upper triangular so that
% QR = Q0*R0 + u*v'.
% GVL4: Section 6.5.1
[m,n] = size(R0);
Q = Q0; R = R0;
w = Q'*u;
% Reduce w to a multiple of e1...
for k=m-1:-1:1
    [c,s] = Givens(w(k),w(k+1));
    J = [c s;-s c];
    idx = k:k+1;
    R(idx,k:n) = J'*R(idx,k:n);
    w(idx)     = J'*w(idx);
    Q(:,idx)   = Q(:,idx)*J;
end
R(1,:) = R(1,:) + w(1)*v';
% Restore R (now upper Hessenberg) to upper triangular form...
for k=1:min(m-1,n)
    [c,s] = Givens(R(k,k),R(k+1,k));
    G = [c s;-s c];
    idx = k:k+1;
    R(idx,k:n) = G'*R(idx,k:n);
    Q(:,idx)   = Q(:,idx)*G;
end
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