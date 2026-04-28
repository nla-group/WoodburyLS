function [Q,R] = QRUpdate2(Q0,R0,U,V,varargin)
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

[m,n] = size(Q0);  r = size(U,2);  U0 = U;
if ~isempty(varargin)
    N = varargin{1};
else
    N = 2;
end

%Orthogonalise U wrt Q0
for i = 1:N  
    U = U-Q0*(Q0'*U);  
end
[Q1,~] = qr(U,0);   
Q = [Q0,Q1];

R = [R0;zeros(r,n)];
W = Q'*U0;
Qv = zeros(1+r,n);

% Reduce W to a multipls of e1,..,er using givens roatation, and apply them
% to Q and R
for i = 1:r
    for k=n-1+r:-1:0+i
        [c,s] = Givens(W(k,i),W(k+1,i));
        J = [c s;-s c];
        idx = [k,k+1];
        R(idx,1:n) = J'*R(idx,1:n);
        W(idx,:) = J'*W(idx,:);
        Q(:,idx)   = Q(:,idx)*J;
    end
end
%R now has r subdiagonals
R= R + W*V'; % adding WV' does not introduce any more fill in

%Now return R to upper triangular form using householder reflectors
%Note that the householder vectors can be stored economically and Q recovered later
for i = 1:n
    [v,beta] = House(R(i:i+r,i));
    Qv(1:r+1,i) = v;
    R(i:i+r,i:n) = R(i:i+r,i:n) - (beta*v)*(v'*R(i:i+r,i:n));
end

Q = Q*ReducedBackAccum(Qv(:,1:n),m);
R = R(1:n,1:n);    Q = Q(1:m,1:n);
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

function Q = ReducedBackAccum(Q_fact,m)
% Explicit formation of an orthogonal matrix from its factored form
% representation. Uses backward accumulation.
[r,n] = size(Q_fact);

if nargin==2 && m>n
    Q = [];
    for j=n:-1:1
        v = [Q_fact(:,j);zeros(n-j,1)];
        beta = 2/(v'*v);
        k = n+r-j-1;
        Q = [1 zeros(1,n-j);zeros(k,1) Q];
        if norm(Q_fact(:,j))>0
            Q = Q - (beta*v)*(v'*Q);
        end
    end
else
    Q = eye(max(n-r,1),max(n-r,1));
    for j=n:-1:1
        v = [Q_fact(:,j);zeros(n-j,1)];
        beta = 2/(v'*v);
        k = n+r-j-1;
        Q = [1 zeros(1,k);zeros(k,1) Q];
        if norm(Q_fact(:,j))>0
            Q = Q - (beta*v)*(v'*Q);
        end
    end
end
end