function plot_sim(samples_sim,hplcparam_sim,metidx,PredType)

tr = samples_sim.(PredType);
y_p_i = prctile(squeeze(tr(:,metidx,:)) ,[5 50 95],1);

tg_i    = hplcparam_sim.tg;
pH_i    = hplcparam_sim.pHo;
mod_i   = hplcparam_sim.Mod;
temp_i  = hplcparam_sim.Temp;
  
subplot(2, 2, 1)
hold on
plotareaprct(pH_i(temp_i==0&mod_i==2&tg_i==30)',y_p_i(:,temp_i==0&mod_i==2&tg_i==30),[0, 0.4470, 0.7410])
plotareaprct(pH_i(temp_i==0&mod_i==2&tg_i==90)',y_p_i(:,temp_i==0&mod_i==2&tg_i==90),[0.6350, 0.0780, 0.1840])
plotareaprct(pH_i(temp_i==0&mod_i==2&tg_i==270)',y_p_i(:,temp_i==0&mod_i==2&tg_i==270),[0.4660, 0.6740, 0.1880])
legend off  
xlim([2 12])

subplot(2, 2, 2)
hold on 
plotareaprct(pH_i(temp_i==0&mod_i==1&tg_i==30)',y_p_i(:,temp_i==0&mod_i==1&tg_i==30),[0, 0.4470, 0.7410])
plotareaprct(pH_i(temp_i==0&mod_i==1&tg_i==90)',y_p_i(:,temp_i==0&mod_i==1&tg_i==90),[0.6350, 0.0780, 0.1840])
plotareaprct(pH_i(temp_i==0&mod_i==1&tg_i==270)',y_p_i(:,temp_i==0&mod_i==1&tg_i==270),[0.4660, 0.6740, 0.1880])   
legend off  
xlim([2 12])

   
subplot(2, 2, 3)
hold on
plotareaprct(pH_i(temp_i==1&mod_i==2&tg_i==30)',y_p_i(:,temp_i==1&mod_i==2&tg_i==30),[0, 0.4470, 0.7410])
plotareaprct(pH_i(temp_i==1&mod_i==2&tg_i==90)',y_p_i(:,temp_i==1&mod_i==2&tg_i==90),[0.6350, 0.0780, 0.1840])
plotareaprct(pH_i(temp_i==1&mod_i==2&tg_i==270)',y_p_i(:,temp_i==1&mod_i==2&tg_i==270),[0.4660, 0.6740, 0.1880])
legend off  
xlim([2 12])

 
subplot(2, 2,4)
hold on 
plotareaprct(pH_i(temp_i==1&mod_i==1&tg_i==30)',y_p_i(:,temp_i==1&mod_i==1&tg_i==30),[0, 0.4470, 0.7410])
plotareaprct(pH_i(temp_i==1&mod_i==1&tg_i==90)',y_p_i(:,temp_i==1&mod_i==1&tg_i==90),[0.6350, 0.0780, 0.1840])
plotareaprct(pH_i(temp_i==1&mod_i==1&tg_i==270)',y_p_i(:,temp_i==1&mod_i==1&tg_i==270),[0.4660, 0.6740, 0.1880])
    
legend off  
xlim([2 12])  

subplot(2,2,1)
hold on
xlim([2 11])
xlabel('pH')
ylabel('t_R')
title(['ACN' ', 25^oC'])
legend off
box off
subplot(2,2,2)
hold on
xlim([2 11])
xlabel('pH')
ylabel('t_R')
title(['MeOH' ', 25^oC'])
legend off
box off
subplot(2,2,3)
hold on
xlim([2 11])
xlabel('pH')
ylabel('t_R')
title(['ACN' ', 35^oC'])
legend off
box off
subplot(2,2,4)
hold on
xlim([2 11])
xlabel('pH')
ylabel('t_R')
title(['MeOH' ', 35^oC'])
legend off
box off
end

function plotareaprct(Time,DVprct,LineColor)

[~,I] = sort(Time);
Time = Time(I);
DVprct = DVprct(:,I);

hold on
fill([Time fliplr(Time)],[DVprct(1,:) fliplr(DVprct(3,:))],LineColor,'LineStyle','none','FaceAlpha',0.4);
plot(Time,DVprct(2,:),'LineWidth',1,'Color',LineColor);
end
