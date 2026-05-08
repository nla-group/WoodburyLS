% demo.m
m = 1e4; n = 500; r = 20;
A = randn(m,n);
U = randn(m,r); V = randn(n,r);
Ahat = A + U*V';
b = Ahat*randn(n ,1);
[Q,R] = qr(A,0);
% solve unmodified LS problem via QR
[x0 , AtAsolver] = WoodburyLS(A,b);
% inefficient : min ||b-(A+UV ')x|| via QR
[Qhat ,Rhat] = qr(Ahat, 0);
x1 = Rhat \(Qhat'*b);
% better to use updated QR factors:
x2 = QRUpdateLS(Q,R,U,V,b);
% even better to solve the modified LS problem like this:
x3 = WoodburyLS(A,b,U,V,x0 , AtAsolver );
% optional : iterative refinement
r3 = b - A*x3 - U*(V'*x3);
x0r = AtAsolver(A'*r3);
e3r = WoodburyLS(A,r3 ,U,V,x0r , AtAsolver );
x3r = x3 + e3r; % refined x3

norm((Ahat*x1-b))  % 1e-11
norm((Ahat*x2-b))  % 1e-11
norm((Ahat*x3-b))  % 1e-9
norm((Ahat*x3r-b)) % 1e-12

[Qhat, Rhat] = QRUpdate(Q,R,U,V);
norm(Qhat*Rhat - Ahat)

[Qhat, Rhat] = Rank1Update(Q,R,U(:,1),V(:,1));
norm(Qhat*Rhat - (A + U(:,1)*V(:,1)'))