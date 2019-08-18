close all;
set(figure,'Position', [50 50 450 200]); hold all; box on;
set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('time'); ylabel('\Delta P');
t_sim=0:0.01:10;
w_bar=3;
u_sim=(sin(t_sim)+sin(3*t_sim+pi/3)+sin(6*t_sim+pi/5)+sin(2*t_sim+pi/5)+(t_sim>1).*(t_sim<3)*100);
u_sim=min(u_sim,w_bar); u_sim=max(u_sim,-w_bar); u_sim=(t_sim>1).*u_sim;
plot(t_sim,u_sim,'b-')
plot(t_sim,w_bar*ones(size(t_sim)),'r--')
plot(t_sim,-w_bar*ones(size(t_sim)),'r--')
ylim([-5 5])

set(figure,'Position', [50 50 450 200]); hold all; box on;
set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('time'); ylabel('\Delta f');
y_sim=lsim(tf([1],[1 1]),u_sim,t_sim);
plot(t_sim,y_sim,'b-')
plot(t_sim,w_bar*ones(size(t_sim)),'r--')
plot(t_sim,-w_bar*ones(size(t_sim)),'r--')
ylim([-5 5])

