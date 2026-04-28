function [Q,R] = QRUpdate3(Q0,R0,U,V)
%Rank-r update of a QR factorisation.
% Uses qr to append orthogonal collumns to Q0.
% Uses qr to reduce W to upper triangular.
% Adds WV^T to R, which causes no additional fill in. 
% Uses qr to reduce R back to upper triangular.
% Q0 is mxn orthogonal, R0 is nxn upper triangular, U is mxr,
% and V is nxr.
% Q is mxn orthogonal and  R is nxn upper triangular so that
% QR = Q0*R0 + U*V'.
n = size(Q0,2);  r = size(U,2);
U1 = U-Q0*(Q0'*U); U2 = U1-Q0*(Q0'*U1);
[Q_,~] = qr(U2,0);   %Orthonormalise U wrt Q0 and save as Q_
Q1 = [Q0,Q_];        %Append Q_ to Q0
[Q2,R2] = qr(Q1'*U); %Calculate Q that upper triangularises Q1'*U
H = Q2'*[R0;zeros(r,n)] + R2*V';  %Form H, weakly uppertrapizoidal 
[Q3,R3] = qr(H);     %Calculate Q that upper triangularises H
Q = Q1*Q2*Q3(:,1:n);  R = R3(1:n,:);%Removes added columns and rows of Q and R.
end