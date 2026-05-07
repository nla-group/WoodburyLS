% demo.m
m = 1e4; n = 500; r = 20;
A = randn(m,n);
U = randn(m,r); V = randn(n,r);
b = (A+U*V')*randn(n ,1) + 1e-6*randn(m,1);
[Q,R] = qr(A);
% solve unmodified LS problem via QR
[x0 , AtAsolver ] = WoodburyLS(A,b);
% inefficient : min ||b-(A+UV ')x|| via QR
Ahat = A + U*V';
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

norm((Ahat*x1-b))
norm((Ahat*x2-b))
norm((Ahat*x3-b))