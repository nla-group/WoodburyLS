##Timings
using LinearAlgebra, BenchmarkTools

include("QRUpdateLS.jl")
include("WoodburyLS.jl")

WoodburyTimes = zeros(10,3);
QRTimes= zeros(10,3);
QRUpdateTimes = zeros(10,3);
m = 100000;
A=[]; b=[]; U=[]; V=[]; solver =[]; x = []; Q = []; R =[]; F = [];  #Must be a better way of doing this..
 for i in 1:10
     for j in 1:1
        n = 100*i;
        r = 100*j
        A = randn(m, n);
        b = randn(m);
        U = randn(m, r);
        V = randn(n, r);
        Q,R = qr(A);
        Q = Matrix(Q);
        # First solve
        # x0 = R\Q'*b;
        # RtR = Cholesky(UpperTriangular(R));
        # AtAsolver = X -> ldiv!(RtR, X);
        x,solver = WoodburyLS_setup(A,b);
        WoodburyTimes[i,j] = @btimed( WoodburyLS(A, b, U, V, x, solver)).time;
        tmp  =@btimed begin
             F = qr(A+U*V'); 
            x = F\b;
        end
        QRTimes[i,j] = tmp.time;
        QRUpdateTimes[i,j] =@btimed( QRUpdateLS(Q,R,U,V,b)).time;
    end
end
 