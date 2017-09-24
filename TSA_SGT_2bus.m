close all; clear;
% Assuming delta_max is between 0 and pi.
M=8; D=2; w=0:0.001:1; Pm=0.2; a=0.8;
delta_eq=asin(Pm/a);
delta_max=pi/2;

figure; hold all; grid on; box on;
delta_plot=-pi:0.01:pi;
plot(delta_plot,sin(delta_plot)-Pm/a)
delta_dash=NR(@(x) sin(x)-sin(delta_eq)-cos(x)*(x-delta_eq),@(x) cos(x)+sin(x)*(x-delta_eq)-cos(x),-pi/2);
delta_min=max(delta_dash,2*delta_eq-pi/2);
g_dash=(sin(delta_min)-sin(delta_eq))/(delta_min-delta_eq);
g_bar=(sin(delta_max)-sin(delta_eq))/(delta_max-delta_eq);
plot(delta_plot,g_dash*(delta_plot-delta_eq),'r--')
plot(delta_plot,g_bar*(delta_plot-delta_eq),'r--')
axis([-pi pi -1.1-Pm/a 1.1-Pm/a])
set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('\delta')

figure;
plot(delta_plot,(sin(delta_plot)-sin(delta_eq))./(delta_plot-delta_eq))


%%
figure; subplot(2,1,1); hold all;
plot(w,sqrt(1./((a*(g_bar+g_dash)/2-M*w.^2).^2+D^2*w.^2)))
scatter(sqrt(a*(g_bar+g_dash)/2/M-D^2/2/M^2),sqrt(1/(D^2*a*(g_bar+g_dash)/2/M-D^4/4/M^2)))
sqrt(1/(D^2*7/18/M-D^4/4/M^2))*11/18;
set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('\omega'); ylabel('\gamma_s');

subplot(2,1,2); hold all;
delta_max=0:0.01:pi;
gamma_s_plot=sqrt(1./(D^2*a*(g_dash+(sin(delta_max)-sin(delta_eq))./(delta_max-delta_eq))/2/M-D^4/4/M^2));
gamma_del_plot=(g_dash-(sin(delta_max)-sin(delta_eq))./(delta_max-delta_eq))/2;
plot(delta_max,gamma_s_plot)
plot(delta_max,gamma_del_plot)
plot(delta_max,gamma_s_plot.*gamma_del_plot)
set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('\delta_{max}'); 
legend('\gamma_s','\gamma_\Delta','\gamma_s\gamma_\Delta'); xlim([0 pi])

%%
figure; hold all;
plot(delta_max,(1-a*gamma_s_plot.*gamma_del_plot)./gamma_s_plot.*delta_max)
set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('\delta_{max}'); ylabel('P_{max}');

% Analytical computation of P_max
g_dash=(sin(delta_min)-sin(delta_eq))/(delta_min-delta_eq);
g_bar=(sin(pi/2-delta_eq)-sin(delta_eq))/(pi/2-delta_eq-delta_eq);
gamma_s=sqrt(1./(D^2*a*(g_dash+g_bar)/2/M-D^4/4/M^2));
gamma_del=(g_dash-g_bar)/2;
P_max=(1-a*gamma_s*gamma_del)/gamma_s*(pi/2-delta_eq);
scatter(pi/2-delta_eq,P_max,100,'r','filled')

% Direct Simulation
del_t_sim=0.1; t_end=200;
t_sim=0:del_t_sim:t_end;
w_dist=sqrt(a*(g_bar+g_dash)/2/M-D^2/2/M^2);
u_intensity=0.1:0.1:0.8; %u_intensity=u_intensity*0.1;
for i=1:size(u_intensity,2)
    x_sim=zeros(2,t_end/del_t_sim);
    x_sim(:,1)=[asin(Pm); 0];
    disturb=u_intensity(i)*sin(w_dist*t_sim);
    
    for t_idx=1:t_end/del_t_sim
        f_sim= @(x) [x(2); 1/M*(-D.*x(2)-a*sin(x(1))+Pm+disturb(t_idx))];
        J_sim=@(x) [0 1; 1/M*(-a*cos(x(1))) -D/M];
        f=@(x) -x+x_sim(:,t_idx)+del_t_sim/2*(f_sim(x_sim(:,t_idx))+f_sim(x));
        J=@(x) -eye(2)+del_t_sim/2*J_sim(x);
        x_sim(:,t_idx+1)=NR(f,J,x_sim(:,t_idx));
    end
    delta_max_sim(i)=max(abs(x_sim(1,:)-delta_eq));
    %figure; hold all; plot(t_sim,x_sim); plot(t_sim,disturb,'r--')
end
plot([0 delta_max_sim],[0 u_intensity]); %xlim([0 pi])

%%
plot_rng=[-1.2*pi 1.2*pi -1 1];
[x1_plot_ij,x2_plot_ij]=meshgrid(linspace(plot_rng(1),plot_rng(2),30),linspace(plot_rng(3),plot_rng(4),30));
for i=1:size(x1_plot_ij,1)
    for j=1:size(x2_plot_ij,2)
        x_plot=[x1_plot_ij(i,j); x2_plot_ij(i,j)];
        dx1_plot(i,j)=x_plot(2);
        dx2_plot(i,j)=1/M*(-D.*x_plot(2)-a*sin(x_plot(1))+Pm);
    end
end
figure; hold all; grid on; box on;
streamslice(x1_plot_ij,x2_plot_ij,dx1_plot,dx2_plot);
plot(x_sim(1,:),x_sim(2,:))

idx_dmax=find(gamma_s_plot.*gamma_del_plot<1);
idx_dmax=idx_dmax(end);
scatter(delta_max(idx_dmax),0,100,'r','filled')

% Future works
% Larger Network
