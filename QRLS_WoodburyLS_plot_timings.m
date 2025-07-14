%% Plot WoodburyLS and QRUpdateLS Timings 



QRTimes = [ 0.137125  0.141676  0.147197;
            0.386224  0.322603  0.371826;
            0.611252  0.51992   0.491516;
            0.731444  0.697242  0.68931;
            1.03788   0.918266  0.920721;
            1.17125   1.13676   1.11172;
            1.42932   1.43361   1.41064;
            1.73168   1.79174   1.91052;
            2.29483   2.44601   2.88365;
            2.77619   3.44304   2.67991];



WoodburyTimes = [0.00465483  0.00648338  0.0102286;
                 0.00582292  0.00864954  0.0136804;
                 0.00741717  0.0108285   0.0169005;
                 0.0094095   0.0139214   0.021239;
                 0.0111224   0.015771    0.0241847;
                 0.0128916   0.0179935   0.0273414;
                 0.014943    0.0214754   0.032032;
                 0.0173255   0.0237776   0.0357641;
                 0.0194048   0.0264571   0.0419747;
                 0.0219504   0.029652    0.0447099];

QRUpdateTimes = [ 0.0178134  0.0324619  0.0570749;
            0.0218296  0.0398709  0.0675537;
            0.028428   0.0472267  0.0765753;
            0.033368   0.0561878  0.088891;
            0.041122   0.066064   0.102405;
            0.0478224  0.077655   0.118085;
            0.0555699  0.0894386  0.136804;
            0.0657897  0.103494   0.155558;
            0.0744778  0.117018   0.178727;
            0.0844221  0.132672   0.199945];


x = linspace(100,1000,10);
figure(1)
clf
mydefaults
hold on 
plot(x,QRTimes(:,1)./WoodburyTimes(:,1),'-x',"Color",[  0    0.4470    0.7410])
plot(x,QRTimes(:,2)./WoodburyTimes(:,2),'-x',"Color",[ 0.8500    0.3250    0.0980])
plot(x,QRTimes(:,3)./WoodburyTimes(:,3),'-x',"Color",[   0.9290    0.6940    0.1250])

plot(x,QRTimes(:,1)./QRUpdateTimes(:,1),'--x',"Color",[  0    0.4470    0.7410])
plot(x,QRTimes(:,2)./QRUpdateTimes(:,2),'--x',"Color",[ 0.8500    0.3250    0.0980])
plot(x,QRTimes(:,3)./QRUpdateTimes(:,3),'--x',"Color",[   0.9290    0.6940    0.1250])
hold off
title('number of rows m = 1e5')
xlabel('number of columns n') 
% legend('r = 10, WoodburyLS','r = 20','r = 30','r = 10, QRUpdate','r = 20','r = 30','Location','NorthWest','FontSize',18)
legend('r = 10','r = 20','r = 30','Location','NorthWest','FontSize',20)
ylabel('Speedup over QR')
grid on
shg
mypdf('WoodburyLS_QRLS_timings_compared',.6,2)

%%
figure(2)
clf
mydefaults
hold on 
plot(x,WoodburyTimes(:,1),'-+',"Color",[  0    0.4470    0.7410])
plot(x,WoodburyTimes(:,2),'-+',"Color",[ 0.8500    0.3250    0.0980])
plot(x,WoodburyTimes(:,3),'-+',"Color",[   0.9290    0.6940    0.1250])

plot(x,QRUpdateTimes(:,1),'-x',"Color",[  0    0.4470    0.7410])
plot(x,QRUpdateTimes(:,2),'-x',"Color",[ 0.8500    0.3250    0.0980])
plot(x,QRUpdateTimes(:,3),'-x',"Color",[   0.9290    0.6940    0.1250])
hold off
title('number of rows m = 1e5')
xlabel('number of columns n') 
legend('r = 10','r = 20','r = 30','Location','NorthWest','FontSize',18)
ylabel('Time in seconds')
grid on
shg
mypdf('WoodburyLS_QRLS_timings',.1,0.8)


