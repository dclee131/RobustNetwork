clear all; close all;
M=8; D=2; w=0:0.001:1; a=0.8; Pm=a*sin(0.2);
delta_eq=asin(Pm/a);
w_max=2;

w_plot=-2*2*pi:0.01:2*2*pi;
delta_cr=NR(@(x) sin(x)-sin(delta_eq)-cos(x)*(x-delta_eq),@(x) cos(x)+sin(x)*(x-delta_eq)-cos(x),-pi/2);
delta_min=max(delta_cr,delta_eq-w_max);
g_underbar=(sin(delta_min)-sin(delta_eq))/(delta_min-delta_eq);
g_bar=(sin(delta_eq+w_max)-sin(delta_eq))/(w_max);

figure; hold all; grid on; box on; set(gcf, 'Position', [100, 100, 400, 300]);
scatter(delta_eq,sin(delta_eq),75,'r','filled');
plot(w_plot,cos(delta_eq)*(w_plot-delta_eq)+Pm/a,'r')
%plot(w_plot,g_underbar*(w_plot-delta_eq)+Pm/a,'b--')
plot(w_plot,(2*cos(delta_eq)-g_bar)*(w_plot-delta_eq)+Pm/a,'b--')
plot(w_plot,g_bar*(w_plot-delta_eq)+Pm/a,'b--')
plot(w_plot,sin(w_plot),'k','LineWidth',1.5)
patch1=patch([delta_eq w_max+delta_eq w_max+delta_eq],[Pm/a g_bar*w_max+Pm/a (2*cos(delta_eq)-g_bar)*w_max+Pm/a],'b','LineStyle','none');
patch2=patch([delta_eq -w_max+delta_eq -w_max+delta_eq],[Pm/a -g_bar*w_max+Pm/a -(2*cos(delta_eq)-g_bar)*w_max+Pm/a],'b','LineStyle','none');
alpha(patch1,0.1); alpha(patch2,0.1)
plot([w_max+delta_eq w_max+delta_eq],[-10 10],'k:','LineWidth',1.8)
plot([-w_max+delta_eq -w_max+delta_eq],[-10 10],'k:','LineWidth',1.8)
x_label=xlabel('$\delta_j-\delta_k$'); set(x_label, 'Interpreter', 'latex','FontSize',15,'FontName','Times New Roman');
y_label=ylabel('$\sin(\delta_j-\delta_k)$'); set(y_label, 'Interpreter', 'latex','FontSize',15,'FontName','Times New Roman');
set(gca,'FontSize',15,'FontName','Times New Roman'); legend('Operating Point','Linearization','Sector Bound','Location','southeast')
axis([-3 3 -1.5 1.5])

%%
figure; hold all; grid on; box on; set(gcf, 'Position', [100, 100, 350, 300]);
global_env=1.255; %1.255
if 0
    patch([0 6 6],[0 -global_env*6 global_env*6],[0 0.4470 0.7410],'LineStyle','--','EdgeColor',[0 0.4470 0.7410],'FaceAlpha',0.2);
    patch([0 -6 -6],[0 -global_env*6 global_env*6],[0 0.4470 0.7410],'LineStyle','--','EdgeColor',[0 0.4470 0.7410],'FaceAlpha',0.2);
end
gain_inf=max(g_underbar-cos(delta_eq),-(g_bar-cos(delta_eq)));
plot(w_plot,sin(w_plot+delta_eq)-Pm/a-cos(delta_eq)*w_plot,'k','LineWidth',2)
if 1
    color_local_sector=[0 0.4470 0.7410];
    plot([w_max w_max],[-10 10],':','color',color_local_sector,'LineWidth',3);
    patch([0 0 0],[0 0 0],1-0.6*(1-color_local_sector),'LineStyle','none');
    %patch0=patch([-10 -10 -10],[-10 -10 -10],1-0.2*(1-[0 0.4470 0.7410]),'LineStyle','--','EdgeColor',[0 0.4470 0.7410]);
    plot([-w_max -w_max],[-10 10],':','color',color_local_sector,'LineWidth',3)
    %patch([0 w_max w_max],[0 -gain_inf*w_max gain_inf*w_max],[0 0.4470 0.7410],'LineStyle','none','FaceAlpha',0.001);
    %patch([0 -w_max -w_max],[0 -gain_inf*w_max gain_inf*w_max],[0 0.4470 0.7410],'LineStyle','none','FaceAlpha',0.001);
    %patch([0 w_max w_max],[0 -gain_inf*w_max gain_inf*w_max],color_local_sector,'LineStyle','none','FaceAlpha',0.3);
    %patch([0 -w_max -w_max],[0 -gain_inf*w_max gain_inf*w_max],color_local_sector,'LineStyle','none','FaceAlpha',0.3);
    patch([0 6 6],[0 -gain_inf*6 gain_inf*6],color_local_sector,'LineStyle','none','FaceAlpha',0.2);
    patch([0 -6 -6],[0 -gain_inf*6 gain_inf*6],color_local_sector,'LineStyle','none','FaceAlpha',0.2);
    %plot(w_plot,-(g_bar-cos(delta_eq))*w_plot,'--','color',color_local_sector,'LineWidth',1.5)
    %plot(w_plot,(g_bar-cos(delta_eq))*w_plot,'--','color',color_local_sector,'LineWidth',1.5)
end
scatter(0,0,75,'r','filled');
set(gca,'FontSize',15,'FontName','Times New Roman');
x_label=xlabel('Transformed state ($z$)'); set(x_label, 'Interpreter', 'latex','FontSize',15,'FontName','Times New Roman');
y_label=ylabel('Nonlinear output ($v$)'); set(y_label, 'Interpreter', 'latex','FontSize',15,'FontName','Times New Roman');
[legend_plot,legend_idx]=legend('$v=\psi(z)$','$\bar{z}$','Sector Bound','Location','northeast');
TextInLegend=findobj(legend_idx, 'type', 'text'); set(TextInLegend, 'Interpreter', 'latex','fontsize',15); set(legend_plot, 'Interpreter', 'latex','fontsize',18,'FontName','Times New Roman');
%axis([-2 2 -0.5 0.5])
axis([-5 5 -7 7])

%% Phase portrait
M=1; D=1.2*0.5; Pm=0.2; a=0.8;
delta_eq=asin(Pm/a);
plot_rng=[-3*pi 3*pi -5 5];
[x1_plot_ij,x2_plot_ij]=meshgrid(linspace(plot_rng(1),plot_rng(2),100),linspace(plot_rng(3),plot_rng(4),100));
dx1_plot_ij=x2_plot_ij;
dx2_plot_ij=M^-1*(-D*x2_plot_ij-a*sin(x1_plot_ij)+Pm);
fx_pre=@(x,u) [x(2); M^-1*(-D*x(2)-a*sin(x(1))+Pm)];
Jx_pre=@(x,u) [0 1; -M^-1*a*sin(x(1)) -M^-1*D*x(2)];
sim_result=zeros(size(size(x1_plot_ij)));
for i=1:size(x1_plot_ij,1)
    for j=1:size(x1_plot_ij,2)
        disp([i,j])
        [x_sim,t_sim]=tds([x1_plot_ij(i,j); x2_plot_ij(i,j)], 0, 5, 0.1, 'trapz', fx_pre, Jx_pre);
        sim_result(i,j)=x_sim(1,end);
    end
end

%%
figure; box on; hold all; set(gcf, 'Position', [100, 100, 450, 350]);
plot(0,0,'color',[0 0.4470 0.7410],'LineWidth',3)
plot(0,0,'r:','LineWidth',3)
streamslice(x1_plot_ij,x2_plot_ij,dx1_plot_ij,dx2_plot_ij);
contour(x1_plot_ij,x2_plot_ij,(sim_result-delta_eq).^2,[1 1]*12,'LineWidth',3,'color',[0 0.4470 0.7410]);
%mesh(x1_plot_ij,x2_plot_ij,(sim_result-delta_eq).^2);
plot([plot_rng(1) plot_rng(2)],pi*[1 1],'r:','LineWidth',3)
plot([plot_rng(1) plot_rng(2)],-pi*[1 1],'r:','LineWidth',3)
scatter(delta_eq,0,'r','filled')
axis(plot_rng);

set(gca,'FontSize',15,'FontName','Times New Roman'); 
xlabel1=xlabel('$\Delta\delta$ (rad)'); set(xlabel1, 'Interpreter', 'latex','FontSize',20);
ylabel1=ylabel('$\dot{\delta}$ (rad/s)'); set(ylabel1, 'Interpreter', 'latex','FontSize',20);
[legend_plot,legend_idx]=legend('Synchronization Constraint','Operational Constraint','Location','northeast');
TextInLegend=findobj(legend_idx, 'type', 'text'); set(TextInLegend, 'Interpreter', 'latex','fontsize',15); set(legend_plot, 'Interpreter', 'latex','fontsize',16,'FontName','Times New Roman');


