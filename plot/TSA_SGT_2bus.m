close all; clear;
M=1; D=1.2; Pm=0.2; a=0.8;
%a=0.1;
delta_eq=asin(Pm/a);

%% Build the model
A=[0 1; -a*cos(delta_eq)/M -D/M];
B_u=[0; 1/M];
B_v=[0; -a/M];
C_y=[0 1];
C_w=[1 0];

G(1,1)=ss(A,B_u,C_y,0);
G(1,2)=ss(A,B_v,C_y,0);
G(2,1)=ss(A,B_u,C_w,0);
G(2,2)=ss(A,B_v,C_w,0);

[y_impulse,t_impulse]=impulse(G(1,1)); gain_mtx(1,1)=trapz(t_impulse,abs(y_impulse));
[y_impulse,t_impulse]=impulse(G(1,2)); gain_mtx(1,2)=trapz(t_impulse,abs(y_impulse));
[y_impulse,t_impulse]=impulse(G(2,1)); gain_mtx(2,1)=trapz(t_impulse,abs(y_impulse));
[y_impulse,t_impulse]=impulse(G(2,2)); gain_mtx(2,2)=trapz(t_impulse,abs(y_impulse));

gamma_del=@(w) cos(delta_eq)-(sin(abs(delta_eq)+w)-sin(abs(delta_eq)))./w;

%% Find Upper bound with simulation
if 1
    % Direct Simulation
    del_t_sim=0.1; t_end=200;
    t_sim=0:del_t_sim:t_end;
    u_intensity=0.1:0.1:1.2;
    %u_intensity=0.5;
    delta_max_sim=zeros(size(u_intensity)); omega_max_sim=zeros(size(u_intensity));
    for i=1:size(u_intensity,2)
        x_sim=zeros(2,t_end/del_t_sim);
        x_sim(:,1)=[delta_eq; 0];
        disturb=u_intensity(i)*ones(size(t_sim));
        
        for t_idx=1:t_end/del_t_sim
            f_sim= @(x) [x(2); 1/M*(-D.*x(2)-a*sin(x(1))+Pm+disturb(t_idx))];
            J_sim=@(x) [0 1; 1/M*(-a*cos(x(1))) -D/M];
            f=@(x) -x+x_sim(:,t_idx)+del_t_sim/2*(f_sim(x_sim(:,t_idx))+f_sim(x));
            J=@(x) -eye(2)+del_t_sim/2*J_sim(x);
            x_sim(:,t_idx+1)=NR(f,J,x_sim(:,t_idx));
        end
        delta_max_sim(i)=max(abs(x_sim(1,:)-delta_eq));
        omega_max_sim(i)=max(abs(x_sim(2,:)));
        %if delta_max_sim(i)>=pi; delta_max_sim(i:end)=1e10; omega_max_sim(i:end)=1e10; break; end
        if size(u_intensity,2)==1; figure; hold all; plot(t_sim,x_sim); plot(t_sim,disturb,'r--'); end
    end
end

%% Worst case input for LTI
if 1
    % Direct Simulation
    del_t_sim=0.1; t_end=200;
    t_sim=0:del_t_sim:t_end;
    [y_impulse,t_impulse]=impulse(G(2,1),t_sim);
    u_intensity=0.1:0.1:1.2;
    %u_intensity=0.5;
    delta_max_sim2=zeros(size(u_intensity)); omega_max_sim2=zeros(size(u_intensity));
    for i=1:size(u_intensity,2)
        x_sim=zeros(2,t_end/del_t_sim);
        x_sim(:,1)=[delta_eq; 0];
        disturb=u_intensity(i)*[ones((size(t_sim)+1)/2) sign(fliplr(y_impulse(1:(size(t_sim,2)-1)/2)'))];
        %disturb=u_intensity(i)*[sign(fliplr(y_impulse(1:(size(t_sim,2)-1)/2)')) ones((size(t_sim)+1)/2)];
        %disturb=u_intensity(i)*sign(fliplr(y_impulse'));
        
        for t_idx=1:t_end/del_t_sim
            f_sim= @(x) [x(2); 1/M*(-D.*x(2)-a*sin(x(1))+Pm+disturb(t_idx))];
            J_sim=@(x) [0 1; 1/M*(-a*cos(x(1))) -D/M];
            f=@(x) -x+x_sim(:,t_idx)+del_t_sim/2*(f_sim(x_sim(:,t_idx))+f_sim(x));
            J=@(x) -eye(2)+del_t_sim/2*J_sim(x);
            x_sim(:,t_idx+1)=NR(f,J,x_sim(:,t_idx));
        end
        delta_max_sim2(i)=max(abs(x_sim(1,:)-delta_eq));
        omega_max_sim2(i)=max(abs(x_sim(2,:)));
        %if delta_max_sim2(i)>=pi; delta_max_sim2(i:end)=1e10; omega_max_sim2(i:end)=1e10; break; end
        if size(u_intensity,2)==1; figure; hold all; plot(t_sim,x_sim); plot(t_sim,disturb,'r--'); end
    end
end

%%
close all;
set(figure,'Position', [50 50 450 450]);
subplot(2,1,1); hold all; box on; grid on; title('(a) Gain with respect to the angle bound')
delta_max_plot=0.01:0.1:pi;
ubar_plot=(1-gain_mtx(2,2).*gamma_del(delta_max_plot)).*delta_max_plot/gain_mtx(2,1);

plot(delta_max_plot,gain_mtx(2,2)*ones(size(delta_max_plot)),'color',[0.9290 0.6940 0.1250],'LineWidth',1.5)
plot(delta_max_plot,gamma_del(delta_max_plot),'color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
plot(delta_max_plot,gain_mtx(2,2).*gamma_del(delta_max_plot),'color',[0 0.4470 0.7410],'LineWidth',1.5)
plot([0 3],[1 1],'k:','LineWidth',2)
plot([2.3873 2.3873],[0 3],'k:','LineWidth',2)
set(gca,'FontSize',15,'FontName','Times New Roman');
xlabel1=xlabel('$\bar{w}$ ($\Delta \delta$ in rad)'); ylabel1=ylabel('$\gamma$');  set(xlabel1, 'Interpreter', 'latex'); set(ylabel1, 'Interpreter', 'latex');
legend1=legend('$\gamma_{w,v}$','$\gamma_\phi$','$\gamma_{w,v}\gamma_\phi$','Orientation','horizontal'); 
set(legend1, 'Interpreter', 'latex','FontSize',17,'FontName','Times New Roman');
axis([0 3 0 1.5]);

subplot(2,1,2); hold all; box on; grid on; title('(b) Perturbation bound')
patch([-1 -1 -1],[-1 -1 -1], 1-0.1*(1-[0 0.4470 0.7410]),'EdgeColor',[0 0.4470 0.7410],'LineWidth',1.5,'facealpha',0.1)
plot([0 delta_max_sim2],[0 u_intensity],'--o','color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
%plot([0 delta_max_sim],[0 u_intensity],'--x','color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
patch([0 delta_max_plot],[0 ubar_plot], [0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410],'LineWidth',1.5,'facealpha',0.1)
plot([0 delta_max_plot],1/gain_mtx(2,1)*[0 delta_max_plot],':','color',[1 0 0]','LineWidth',1.5)
plot([2.3873 2.3873],[0 3],'k:','LineWidth',2)
scatter(1.21,0.501,50,'r','filled')
set(gca,'FontSize',15,'FontName','Times New Roman'); 
xlabel1=xlabel('$\bar{w}$ ($\Delta \delta$ in rad)'); ylabel1=ylabel('$\bar{u}$ ($\Delta P$ in p.u.)'); set(xlabel1, 'Interpreter', 'latex'); set(ylabel1, 'Interpreter', 'latex');
legend1=legend('$\gamma_{w,u}\bar{u}\leq(I-\gamma_{w,v}\gamma_\varphi(\bar{w}))\bar{w}$','Simulations','Orientation','horizontal');
set(legend1, 'Interpreter', 'latex','FontSize',15,'FontName','Times New Roman');
axis([0 3 0 0.8]);

set(figure,'Position', [50 50 450 250]); hold all; box on; grid on;
ybar_plot=gain_mtx(1,1)*ubar_plot+gain_mtx(1,2)*gamma_del(delta_max_plot).*delta_max_plot;
patch([-1 -1 -1],[-1 -1 -1], 1-0.1*(1-[0 0.4470 0.7410]),'EdgeColor',[0 0.4470 0.7410],'LineWidth',1.5,'facealpha',0.1)
patch([0 u_intensity 10],[0 omega_max_sim2 -1]/(2*pi), [0.8500 0.3250 0.0980],'Marker','o','LineStyle','--','EdgeColor',[0.8500 0.3250 0.0980],'LineWidth',1.5,'facealpha',0.1)
%plot([0 u_intensity],[0 omega_max_sim2]/(2*pi),'--o','color',[0.8500 0.3250 0.0980]','LineWidth',1.5)
%%plot([0 u_intensity],[0 omega_max_sim]/(2*pi),'--x','color',[0.8500 0.3250 0.0980]','LineWidth',1.5)
%patch([0 u_intensity 10],gain_mtx(1,1)*[0 u_intensity -1]/(2*pi), [1 0 0],'LineStyle','--','EdgeColor',[1 0 0],'LineWidth',1.5,'facealpha',0.1)
plot([0 u_intensity],gain_mtx(1,1)*[0 u_intensity]/(2*pi),':','color',[1 0 0]','LineWidth',2)
%patch([0 ubar_plot],[0 ybar_plot]/(2*pi), [0 0.4470 0.7410],'EdgeColor',[0 0.4470 0.7410],'LineWidth',1.5,'facealpha',0.1)
%%plot([0 5],[0.1 0.1],'k:','LineWidth',2); scatter(0.435,0.1,50,'r','filled')
set(gca,'FontSize',15,'FontName','Times New Roman'); 
xlabel1=xlabel('Disturbance, $\bar{u}$ ($\Delta P$ in p.u.)'); set(xlabel1, 'Interpreter', 'latex');
ylabel1=ylabel('Output, $\bar{y}$ ($\dot{\delta}$ in Hz)'); set(ylabel1, 'Interpreter', 'latex');
%legend1=legend('$\gamma_{y,u}\bar{u}+\gamma_{y,v}\gamma_\varphi(\bar{w})\bar{w}\leq\bar{y}$','Simulations','Linearization','Orientation','vertical');
set(legend1, 'Interpreter', 'latex','FontSize',15,'FontName','Times New Roman');
axis([0 0.8 0 0.15]);

return

%%
%sys=tf([a 0],[M D a*cos(delta_eq)]); [y,t]=impulse(sys,100);
%gamma_freq_plot=trapz(t,abs(y));
%delta_max=delta_max_plot(idx_delmax);
%freq_max=gamma_freq_plot*gamma_del_plot(idx_delmax)*delta_max;

plot_rng=[-1.2*pi 1.2*pi -1 1];
[x1_plot_ij,x2_plot_ij]=meshgrid(linspace(plot_rng(1),plot_rng(2),100),linspace(plot_rng(3),plot_rng(4),100));
for i=1:size(x1_plot_ij,1)
    for j=1:size(x2_plot_ij,2)
        x_plot=[x1_plot_ij(i,j); x2_plot_ij(i,j)];
        dx1_plot(i,j)=x_plot(2);
        dx2_plot(i,j)=1/M*(-D*x_plot(2)-a*sin(x_plot(1))+Pm);
    end
end

if 0
    figure; hold all; grid on; box on;
    streamslice(x1_plot_ij,x2_plot_ij,dx1_plot,dx2_plot);
    patch([delta_eq+delta_max delta_eq+delta_max delta_eq-delta_max delta_eq-delta_max],[freq_max -freq_max -freq_max freq_max],'blue')
    alpha(0.1)
    axis(plot_rng)
    set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('\delta'); ylabel('\omega');
end
return;

%% Region of Attraction
figure; box on; hold all;
sys=tf(M,[M D a*cos(delta_eq)]); [y_impulse,t_impulse]=impulse(sys,300);
freq_gain=max(abs(y_impulse));

for i=1:size(x1_plot_ij,1)
    for j=1:size(x2_plot_ij,2)
        x_plot=[x1_plot_ij(i,j); x2_plot_ij(i,j)];
        delta_ij=abs(x_plot(1)-delta_eq); freq_ij=abs(x_plot(2));
        delta_rng=linspace(0,pi/2,100);
        u_max=max(delta_rng-gain_mtx(2,2)*gamma_del(delta_rng).*delta_rng);
        feas_region(i,j)=u_max-delta_ij-freq_gain*freq_ij;
    end
end

streamslice(x1_plot_ij,x2_plot_ij,dx1_plot,dx2_plot);
contour(x1_plot_ij,x2_plot_ij,feas_region,[0 0],'LineWidth',3)
%mesh(x1_plot_ij,x2_plot_ij,feas_region)
axis(plot_rng)