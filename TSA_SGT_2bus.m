close all; clear;
% Assuming delta_max is between 0 and pi.
M=8; D=2; w=0:0.001:1; Pm=0.2; a=0.8;
delta_eq=asin(Pm/a);
delta_max=pi/2;

delta_dash=NR(@(x) sin(x)-sin(delta_eq)-cos(x)*(x-delta_eq),@(x) cos(x)+sin(x)*(x-delta_eq)-cos(x),-pi/2);
delta_min=max(delta_dash,2*delta_eq-pi/2);
g_dash=(sin(delta_min)-sin(delta_eq))/(delta_min-delta_eq);
g_bar=(sin(delta_max)-sin(delta_eq))/(delta_max-delta_eq);

%%
subplot(2,1,1); hold all; box on; grid on;
delta_max_plot=0.01:0.01:2;
sys=tf(1,[M D a*cos(delta_eq)]); [y_impulse,t_impulse]=impulse(sys,300);
gamma_s_plot=trapz(t_impulse,abs(y_impulse))*ones(size(delta_max_plot));
g_bar=abs((sin(delta_max_plot+delta_eq)-sin(delta_eq))./delta_max_plot-cos(delta_eq));
gamma_del_plot=g_bar;
%g_dash=abs((sin(delta_eq)-sin(delta_max_plot-delta_eq))./(delta_max_plot)-cos(delta_eq));
%gamma_del_plot=max(g_bar,g_dash);
plot(delta_max_plot,gamma_s_plot)
plot(delta_max_plot,gamma_del_plot)
plot(delta_max_plot,a*gamma_s_plot.*gamma_del_plot)
set(gca,'FontSize',15,'FontName','Times New Roman'); ylabel('\gamma')
legend('\gamma_G','\gamma_\Delta','\gamma_s\gamma_\Delta','Orientation','horizontal'); 
axis([0 2 0 3])

subplot(2,1,2); hold all; box on; grid on;
plot(delta_max_plot,(1-a*gamma_s_plot.*gamma_del_plot)./gamma_s_plot.*delta_max_plot)
set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('|\delta_{max}-\delta_0|'); ylabel('\Delta P_{max}');
axis([0 2 0 0.65])

%% Upper bound
if 0
    % Direct Simulation
    del_t_sim=0.1; t_end=200;
    t_sim=0:del_t_sim:t_end;
    u_intensity=0.1:0.1:0.8; %u_intensity=u_intensity*0.1;
    for i=1:size(u_intensity,2)
        x_sim=zeros(2,t_end/del_t_sim);
        x_sim(:,1)=[asin(Pm); 0];
        %disturb=u_intensity(i)*sin(w_dist*t_sim);
        disturb=u_intensity(i)*ones(size(t_sim));
        
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
    plot([0 delta_max_sim+1e10*(delta_max_sim>pi)],[0 u_intensity]); %xlim([0 pi])
    legend('Lower bound','Upper bound','Orientation','horizontal')
end

%%
del_t_sim=0.1; t_end=300; t_sim=0:del_t_sim:t_end;
[Pmax, idx_delmax]=max((1-a*gamma_s_plot.*gamma_del_plot)./gamma_s_plot.*delta_max_plot);
for idx_case=1
    x_sim{idx_case}=zeros(2,t_end/del_t_sim);
    x_sim{idx_case}(:,1)=[asin(Pm); 0];
    %if idx_case==1; disturb=0.5*ones(size(t_sim)); end
    %if idx_case==2; disturb=-0.7*ones(size(t_sim)); end
    %if idx_case==3; disturb=1*sign(interp1(t_impulse,y_impulse,t_sim)); end
    %if idx_case==4; disturb=-0.5*sign(interp1(t_impulse,y_impulse,t_sim)); end
    for t_idx=1:t_end/del_t_sim
        disturb(t_idx)=0.6*sign(x_sim{idx_case}(1,t_idx)-delta_eq);
        f_sim= @(x) [x(2); 1/M*(-D.*x(2)-a*sin(x(1))+Pm+disturb(t_idx))];
        J_sim=@(x) [0 1; 1/M*(-a*cos(x(1))) -D/M];
        f=@(x) -x+x_sim{idx_case}(:,t_idx)+del_t_sim/2*(f_sim(x_sim{idx_case}(:,t_idx))+f_sim(x));
        J=@(x) -eye(2)+del_t_sim/2*J_sim(x);
        x_sim{idx_case}(:,t_idx+1)=NR(f,J,x_sim{idx_case}(:,t_idx));
    end
end

sys=tf([a 0],[M D a*cos(delta_eq)]); [y,t]=impulse(sys,100);
gamma_freq_plot=trapz(t,abs(y));
delta_max=delta_max_plot(idx_delmax);
freq_max=gamma_freq_plot*gamma_del_plot(idx_delmax)*delta_max;

plot_rng=[-1.2*pi 1.2*pi -1 1];
[x1_plot_ij,x2_plot_ij]=meshgrid(linspace(plot_rng(1),plot_rng(2),30),linspace(plot_rng(3),plot_rng(4),30));
for i=1:size(x1_plot_ij,1)
    for j=1:size(x2_plot_ij,2)
        x_plot=[x1_plot_ij(i,j); x2_plot_ij(i,j)];
        dx1_plot(i,j)=x_plot(2);
        dx2_plot(i,j)=1/M*(-D*x_plot(2)-a*sin(x_plot(1))+Pm);
    end
end

if 1
    figure; hold all; grid on; box on;
    streamslice(x1_plot_ij,x2_plot_ij,dx1_plot,dx2_plot);
    plot(x_sim{1}(1,:),x_sim{1}(2,:),'b','LineWidth',2);
    %plot(x_sim{2}(1,:),x_sim{2}(2,:),'r','LineWidth',2)
    %plot(x_sim{3}(1,:),x_sim{1}(2,:),'g','LineWidth',2); plot(x_sim{4}(1,:),x_sim{2}(2,:),'m','LineWidth',2)
    patch([delta_eq+delta_max delta_eq+delta_max delta_eq-delta_max delta_eq-delta_max],[freq_max -freq_max -freq_max freq_max],'blue')
    alpha(0.1)
    axis(plot_rng)
    set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('\delta'); ylabel('\omega');
end