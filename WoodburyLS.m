function [x,AtAsolver] = WoodburyLS(A,b,U,V,x0,AtAsolver)
%WoodburyLS Solves the least squares problem min_x ||b - (A+UV')x||
% where A is an m x n matrix, m >= n, and U and V have r columns.
% Requires A and A+UV' to be full rank.
%
% First call:
%   [x0,AtAsolver] = WoodburyLS(A,b) returns LS solution x0 = A\b
%   and a function AtAsolver that solves A'A x = b for a given b.
%
% Every follow-up call:
%   x = WoodburyLS(A,b,U,V,x0,AtAsolver) returns the LS solution
%   x = (A+U*V')\b using a given solution x0 = A\b and AtAsolver.

if nargin < 3
    [Q,R] = qr(A,0); x = R \ (Q'*b);
    AtAsolver = @(X) R \ (R' \ X);
    return
end
r = size(U,2);
X  = [V, A'*U];                       % X = [V, A'U]
Yt = [X(:,r+1:2*r)' + (U'*U)*V'; V']; % Y = [(A+UV')U, V]
Z = AtAsolver(X);                     % Z = (A'A)\X
w = x0 + Z(:, 1:r) * (U'*b);          % Z(:,1:r) = (A'A)\V
Mw = Z * ((eye(2*r)+Yt*Z) \ (Yt*w));  % M = (A'A)\X(I+Y'(A'A)\X)\Y'
x = w - Mw;