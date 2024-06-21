close all
clear all
clc

runs = 10;

m = 1e5; 
NN = [ 100:100:1000 ];
RR = [10,20,30];

T_qr0 = [];
T_qr1 = [];
T_chol = [];
T_update = [];
ACC = [];

for i = 1:length(NN), n = NN(i); 
    disp(['cols n = ' num2str(n) ])
    for j = 1:length(RR), r = RR(j);
        disp(['rank n = ' num2str(r) ])
        rng('default')
        A = randn(m,n); b = randn(m,1);
        disp('solve min || b - Ax || using QR')
        
        tt = [];
        for run = 1:runs
            tic
            [x0,AtAsolver] = WoodburyLS(A,b); % original LS problem        
            tt(end+1) = toc;
        end
        T_qr0(i,j,:) = tt;

        U = randn(m,r); V = randn(n,r);
        disp('solve min || b - (A+UV'')x || using QR')
        tt = [];
        for run = 1:runs
            tic
            Ahat = A + U*V';
            [Qhat,Rhat] = qr(Ahat,0);
            x1 = Rhat\(Qhat'*b);    % updated problem (from scratch)
            tt(end+1) = toc;
        end
        T_qr1(i,j,:) = tt;

        disp('solve min || b - (A+UV'')x || using Chol-QR')
        tt = [];
        for run = 1:runs
            tic
            Ahat = A + U*V';
            M = Ahat'*Ahat;
            Rc = chol(M);
            %Qc = Ahat/Rc; %xc = Rc\(Qc'*b);
            xc = Rc\((Rc')\(Ahat'*b));
            tt(end+1) = toc;
        end
        T_chol(i,j,:) = tt;
        
        % compute update like this instead:
        disp('solve min || b - (A+UV'')x || using update formula')
        tt = [];
        for run = 1:runs
            tic
            x2 = WoodburyLS(A,b,U,V,x0,AtAsolver);
            tt(end+1) = toc;
        end
        T_update(i,j,:) = tt;

        ERR(i,j) = norm(x2 - x1)/norm(x1);      % error check 
        
        save timings NN RR T_qr0 T_qr1 T_chol T_update ERR    
    end
end


%%
mydefaults
load timings
plot(NN,mean(T_qr1,3)./mean(T_update,3),'-+')
title('number of rows m = 1e5')
xlabel('number of columns n') 
legend('r = 10','r = 20','r = 30','Location','NorthWest','FontSize',18)
ylabel('speedup over QR')
grid on
shg
mypdf('WoodburyLS_timings',.6,0.8)


