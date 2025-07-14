##Timings
using LinearAlgebra, BenchmarkTools

include("QRUpdateLS.jl")
include("WoodburyLS.jl")

Ns = 10; Rs = 3;
WoodburyTimes = zeros(Ns,Rs);
QRTimes= zeros(Ns,Rs);
QRUpdateTimes = zeros(Ns,Rs);
m = 100000;
A=[]; b=[]; U=[]; V=[]; solver =[]; x = []; Q = []; R =[]; F = [];  #Must be a better way of doing this..
erW = zeros(Ns,Rs);
erQR = zeros(Ns,Rs);
 for i in 1:Ns
     for j in 1:Rs
        n = 100*i;  r = 10*j
        A = randn(m, n); b = randn(m);
        U = randn(m, r); V = randn(n, r);
        Q,R = qr(A);  Q = Matrix(Q);
        # First solve
        x,solver = WoodburyLS_setup(A,b);
        tmp = @btimed( WoodburyLS(A, b, U, V, x, solver));
        WoodburyTimes[i,j] =tmp.time;
        x1 = tmp.value[1];
        tmp  =@btimed begin
             F = qr(A+U*V'); 
            x0 = F\b;
        end
        x0 = tmp.value[1];
        QRTimes[i,j] = tmp.time;
        tmp = @btimed( QRUpdateLS(Q,R,U,V,b));
        QRUpdateTimes[i,j] = tmp.time;
        x2 = tmp.value[1]
        # x1 = WoodburyLS(A, b, U, V, x, solver)
        erW[i,j] = norm(x1-x0)/norm(x0)
        # x2  = QRUpdateLS(Q,R,U,V,b)
        erQR[i,j] = norm(x2-x0)/norm(x0)
    end
end
 

# Iterative refinement
tmp = @btimed begin
r1 = b - A*x1 - U*(V'*x1);                     #   
xr = solver(A'*r1);                            #0.0383605 seconds
e1r = WoodburyLS(A, r1, U, V, xr, solver);     #
x1r = x1 + e1r;                                # refined x1       
end