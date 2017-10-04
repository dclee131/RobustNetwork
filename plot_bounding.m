clear all; close all;
M=8; D=2; w=0:0.001:1; a=0.8; Pm=a*sin(1);
delta_eq=asin(Pm/a);
w_max=1;

w_plot=-2*pi:0.01:2*pi;
delta_cr=NR(@(x) sin(x)-sin(delta_eq)-cos(x)*(x-delta_eq),@(x) cos(x)+sin(x)*(x-delta_eq)-cos(x),-pi/2);
delta_min=max(delta_cr,delta_eq-w_max);
g_underbar=(sin(delta_min)-sin(delta_eq))/(delta_min-delta_eq);
g_bar=(sin(delta_eq+w_max)-sin(delta_eq))/(w_max);

figure; hold all; grid on; box on; set(gcf, 'Position', [100, 100, 400, 300]);
scatter(delta_eq,sin(delta_eq),75,'r','filled');
plot(w_plot,cos(delta_eq)*(w_plot-delta_eq)+Pm/a,'r')
plot(w_plot,g_underbar*(w_plot-delta_eq)+Pm/a,'b--')
plot(w_plot,g_bar*(w_plot-delta_eq)+Pm/a,'b--')
plot(w_plot,sin(w_plot),'k','LineWidth',1)
%patch1=patch([delta_eq w_max+delta_eq w_max+delta_eq],[Pm/a g_bar*w_max+Pm/a g_underbar*w_max+Pm/a],'b','LineStyle','none');
%patch2=patch([delta_eq -w_max+delta_eq -w_max+delta_eq],[Pm/a -g_bar*w_max+Pm/a -g_underbar*w_max+Pm/a],'b','LineStyle','none');
%alpha(patch1,0.3); alpha(patch2,0.3)
plot([w_max+delta_eq w_max+delta_eq],[-10 10],'k:','LineWidth',1.8)
plot([-w_max+delta_eq -w_max+delta_eq],[-10 10],'k:','LineWidth',1.8)
set(gca,'FontSize',15,'FontName','Times New Roman'); legend('Nominal Operating Point','Linearization','Sector Bound','Location','southeast')
axis([-3 3 -1.5 1.5])


%%
figure; hold all; grid on; box on; set(gcf, 'Position', [100, 100, 400, 300]);
scatter(0,0,75,'r','filled');
plot(w_plot,sin(w_plot+delta_eq)-Pm/a-cos(delta_eq)*w_plot,'k','LineWidth',1)
plot(w_plot,(g_underbar-cos(delta_eq))*w_plot,'b--')
plot(w_plot,(g_bar-cos(delta_eq))*w_plot,'b--')
gain_inf=max(g_underbar-cos(delta_eq),-(g_bar-cos(delta_eq)));
patch1=patch([0 w_max w_max],[0 -gain_inf*w_max gain_inf*w_max],'b','LineStyle','none');
patch2=patch([0 -w_max -w_max],[0 -gain_inf*w_max gain_inf*w_max],'b','LineStyle','none');
alpha(patch1,0.2); alpha(patch2,0.2)
plot([w_max w_max],[-10 10],'k:','LineWidth',1.8)
plot([-w_max -w_max],[-10 10],'k:','LineWidth',1.8)
set(gca,'FontSize',15,'FontName','Times New Roman'); %xlabel('w'); ylabel('v');
axis([-2 2 -1 1])

