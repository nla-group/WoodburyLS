function [Q,R] = Rank1Update(Q0,R0,u,v)
% Rank-1 update of a QR factorization.
% Q0 is mxn orthogonal, R0 is mxn upper triangular, u is mx1,
% and v is nx1.
% Q is mxm orthogonal and R is mxn upper triangular so that
% QR = Q0*R0 + u*v'.
% This is a slightly modified version of
% https://www.cs.cornell.edu/cv/GVL4/M-Files/Chapter%206/Functions/Rank1Update.m
% as described in Golub--Van Loan 4th Ed: Section 6.5.1
[m,n] = size(R0);
Q = Q0; R = R0;
w = Q'*u;
% Reduce w to a multiple of e1...
for k=m-1:-1:1
    idx = [k,k+1];
    G = planerot(w(idx));
    R(idx,:) = G*R(idx,:);
    w(idx)   = G*w(idx);
    Q(:,idx) = Q(:,idx)*G';
end
R(1,:) = R(1,:) + w(1)*v';
% Restore R (now upper Hessenberg) to upper triangular form...
for k=1:min(m-1,n)
    idx = [k,k+1];
    G = planerot(R(idx,k));
    R(idx,k:n) = G*R(idx,k:n);
    Q(:,idx)   = Q(:,idx)*G';
end
end