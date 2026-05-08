function [Q,R] = QRUpdate(Q0,R0,U,V,varargin)
%Rank-r update of a QR factorisation.
% Uses qr to append orthogonal collumns to Q0.
% Applys Givens Givens rotations to reduce W to upper triangular,
% introducing r diagonals below the main diagonal of R.
% Adds WV^T to R, which causes no additional fill in.
% Uses householder reflections to reduce R back to upper triangular.
% Q0 is mxn orthogonal, R0 is nxn upper triangular, U is mxr,
% and V is nxr.
% Q is mxn orthogonal and  R is nxn upper triangular so that
% QR = Q0*R0 + U*V'.
% Bloor Riley et al, A Sherman–Morrison–Woodbury Approach to Solving
% Least Squares Problems with Low-Rank Updates, 2026.

[m,n] = size(Q0);  r = size(U,2);
% orthogonalise U wrt Q0
U1 = U-Q0*(Q0'*U);
U1 = U1-Q0*(Q0'*U1); % perform orthogonalisation twice

[Q1,~] = qr(U1,0);
Q = [Q0,Q1];
R = [R0;zeros(r,n)];
W = Q'*U0;

% reduce W to upper triangular using Givens rotations, and apply to R and Q
for i = 1:r
    for k=m-1:-1:i
        idx = [k,k+1];
        G = planerot(W(idx,i));
        R(idx,:) = G*R(idx,:);
        W(idx,:) = G*W(idx,:);
        Q(:,idx) = Q(:,idx)*G';
    end
end
R = R + W*V'; % R is now upper trapezoidal
% now return R to upper triangular form using householder reflectors
for i = 1:n
    [v,beta] = gallery("house", R(i:i+r,i));
    Q(:,i:i+r) = Q(:,i:i+r) - beta*(Q(:,i:i+r)*v)*v';
    R(i:i+r,i:n) = R(i:i+r,i:n) - (beta*v)*(v'*R(i:i+r,i:n));
end
R = R(1:n,1:n);
Q = Q(1:m,1:n);
end