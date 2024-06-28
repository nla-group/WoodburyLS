# WoodburyLS
**Fast solution of least squares problems with low-rank modification**

The MATLAB function `WoodburyLS` allows the fast solution of a least squares problem that underwent a low-rank update. The below code demonstrates its use.

```Matlab
m = 1e5; n = 500; r = 20;
A = randn(m,n); b = randn(m,1);
U = randn(m,r); V = randn(n,r);

% solve unmodified LS problem via QR
[x0,AtAsolver] = WoodburyLS(A,b);              % 2.340 seconds

% inefficient: min ||b-(A+UV')x|| via QR
Ahat = A + U*V';                               %  
[Qhat,Rhat] = qr(Ahat,0);                      % 2.390 seconds
x1 = Rhat\(Qhat'*b);                           % 

% better solve modified LS problem like this:
x2 = WoodburyLS(A,b,U,V,x0,AtAsolver);         % 0.037 seconds
```

**Preprint** 

For more details, check out https://arxiv.org/abs/2406.15120

**BibTeX:**
```
@techreport{GNWB24,
  title   = {A {S}herman--{M}orrison--{W}oodbury approach to solving least squares problems with low-rank updates},
  author  = {G\"{u}ttel, Stefan and Nakatsukasa, Yuji and Webb, Marcus and Bloor Riley, Alban},
  year    = {2024},
  number  = {arXiv:2406.15120},
  pages   = {6},
  institution = {arXiv}, 
  type    = {arXiv EPrint},
  url     = {https://arxiv.org/abs/2406.15120}
}
```
