function plot_data(data,metidx)

data = sortrows(data,{'pHs'},{'ascend'});

ktore1 = (data.METID==metidx)&(data.Temp==25)&(data.Mod=='ACN');
ktore2 = (data.METID==metidx)&(data.Temp==25)&(data.Mod=='MeOH');
ktore3 = (data.METID==metidx)&(data.Temp==35)&(data.Mod=='ACN');
ktore4 = (data.METID==metidx)&(data.Temp==35)&(data.Mod=='MeOH');

subplot(2,2,1)
hold on
plot(data.pHo(ktore1&data.tg==30),data.RT(ktore1&data.tg==30),'.','MarkerSize',15,'Color',[0, 0.4470, 0.7410]);
plot(data.pHo(ktore1&data.tg==90),data.RT(ktore1&data.tg==90),'.','MarkerSize',15,'Color',[0.6350, 0.0780, 0.1840]);
plot(data.pHo(ktore1&data.tg==270),data.RT(ktore1&data.tg==270),'.','MarkerSize',15,'Color',[0.4660, 0.6740, 0.1880]);
subplot(2,2,2)
hold on
plot(data.pHo(ktore2&data.tg==30),data.RT(ktore2&data.tg==30),'.','MarkerSize',15,'Color',[0, 0.4470, 0.7410]);
plot(data.pHo(ktore2&data.tg==90),data.RT(ktore2&data.tg==90),'.','MarkerSize',15,'Color',[0.6350, 0.0780, 0.1840]);
plot(data.pHo(ktore2&data.tg==270),data.RT(ktore2&data.tg==270),'.','MarkerSize',15,'Color',[0.4660, 0.6740, 0.1880]);
subplot(2,2,3)
hold on
plot(data.pHo(ktore3&data.tg==30),data.RT(ktore3&data.tg==30),'.','MarkerSize',15,'Color',[0, 0.4470, 0.7410]);
plot(data.pHo(ktore3&data.tg==90),data.RT(ktore3&data.tg==90),'.','MarkerSize',15,'Color',[0.6350, 0.0780, 0.1840]);
plot(data.pHo(ktore3&data.tg==270),data.RT(ktore3&data.tg==270),'.','MarkerSize',15,'Color',[0.4660, 0.6740, 0.1880]);
subplot(2,2,4)
hold on
plot(data.pHo(ktore4&data.tg==30),data.RT(ktore4&data.tg==30),'.','MarkerSize',15,'Color',[0, 0.4470, 0.7410]);
plot(data.pHo(ktore4&data.tg==90),data.RT(ktore4&data.tg==90),'.','MarkerSize',15,'Color',[0.6350, 0.0780, 0.1840]);
plot(data.pHo(ktore4&data.tg==270),data.RT(ktore4&data.tg==270),'.','MarkerSize',15,'Color',[0.4660, 0.6740, 0.1880]);


subplot(2,2,1)
hold on
xlabel('pH')
ylabel('t_R')
title(['ACN' ', 25^oC'])
xlim([2 11])
legend off
box off
subplot(2,2,2)
hold on
xlabel('pH')
ylabel('t_R')
title(['MeOH' ', 25^oC'])
xlim([2 11])
legend off
box off
subplot(2,2,3)
hold on
xlabel('pH')
ylabel('t_R')
title(['ACN' ', 35^oC'])
legend off
box off
xlim([2 11])
subplot(2,2,4)
xlabel('pH')
ylabel('t_R')
title(['MeOH' ', 35^oC'])
legend off
box off
xlim([2 11])
end