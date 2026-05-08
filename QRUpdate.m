function [Q,R] = QRUpdate(Q0,R0,U,V)
%Rank-r update of a QR factorisation.
% Uses qr to append orthogonal collumns to Q0.
% Applys Givens Givens rotations to reduce W to upper triangular,
% introducing r diagonals below the main diagonal of R.
% Adds WV^T to R, which causes no additional fill in.
% Uses householder reflections to reduce R back to upper triangular.
% Q0 is mxn orthogonal, R0 is nxn upper triangular, U is mxr,
% and V is nxr. We assume m >= n+r.
% Q is mxn orthogonal and  R is nxn upper triangular so that
% QR = Q0*R0 + U*V'.
% Bloor Riley et al, A Sherman–Morrison–Woodbury Approach to Solving
% Least Squares Problems with Low-Rank Updates, 2026.

[n,r] = size(V);
% orthogonalise U wrt Q0
U1 = U-Q0*(Q0'*U);
Qappend = orth(U1,10*eps);
Q = [Q0,Qappend];
R = [R0;zeros(size(Qappend,2),n)];
W = Q'*U;
% reduce W to upper triangular using Givens rotations, and apply to R and Q
for i = 1:r
    for k=size(W,1)-1:-1:i
        idx = [k,k+1];
        G = planerot(W(idx,i));
        R(idx,:) = G*R(idx,:);
        W(idx,:) = G*W(idx,:);
        Q(:,idx) = Q(:,idx)*G';
    end
end
R = R + W*V'; % R is now upper trapezoidal
% now return R to upper triangular form using householder reflectors
for k = 1:(size(R,1)-r)
    [v,beta] = gallery("house", R(k:k+r,k));
    R(k:k+r,k:n) = R(k:k+r,k:n) - (beta*v)*(v'*R(k:k+r,k:n));
    Q(:,k:k+r) = Q(:,k:k+r) - beta*(Q(:,k:k+r)*v)*v';
end
R = R(1:n,1:n);
Q = Q(:,1:n);
end