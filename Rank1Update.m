function [Q,R] = Rank1Update(Q0,R0,u,v)
% Rank-1 update of a QR factorization.
% Q0 is mxn orthogonal, R0 is nxn upper triangular, u is mx1,
% and v is nx1, with m > n.
% Q is mxn orthogonal and R is mxn upper triangular so that
% QR = Q0*R0 + u*v'.
% This is a modified version of
% https://www.cs.cornell.edu/cv/GVL4/M-Files/Chapter%206/Functions/Rank1Update.m
% which requires square Q0.
n = size(Q0,2);
% orthogonalise u wrt Q0
u1 = u-Q0*(Q0'*u);
if norm(u1) > 10*eps
    Q = [Q0,u1/norm(u1)];
    R = [R0;zeros(1,n)];
else
    Q = Q0; R = R0;
end
w = Q'*u;
% Reduce w to a multiple of e1...
for k=length(w)-1:-1:1
    idx = [k,k+1];
    G = planerot(w(idx));
    R(idx,:) = G*R(idx,:);
    w(idx)   = G*w(idx);
    Q(:,idx) = Q(:,idx)*G';
end
R(1,:) = R(1,:) + w(1)*v';
% Restore R (now upper Hessenberg) to upper triangular form...
for k=1:size(R,1)-1
    idx = [k,k+1];
    G = planerot(R(idx,k));
    R(idx,k:n) = G*R(idx,k:n);
    Q(:,idx)   = Q(:,idx)*G';
end
Q = Q(:,1:n);
R = R(1:n,1:n);
end