clear; close all;

sys_case=39;
model_order='2gov'; activate_Lossless=1;
run(['dyn' int2str(sys_case)]) % Get info from data file
run(['model_' model_order]) % Build model

del_t=0.01; t_end=30;
t_sim=0:del_t:t_end;
disturb_mag=0.939;
disturb_bus=[3 15 27];
disturb_cont=zeros(num_var,size(t_sim,2));
disturb_step=zeros(num_var,size(t_sim,2));
for i=1:size(disturb_bus,2); 
    rand_sig=sum(cos(2*rand(3,1)*t_sim+2*pi*rand(3,1))); 
    disturb_cont(disturb_bus(i)+2*num_gen,:)=disturb_mag*rand_sig/max(abs(rand_sig));
    disturb_step(disturb_bus(i)+2*num_gen,t_sim<15)=-disturb_mag;
    disturb_step(disturb_bus(i)+2*num_gen,t_sim>15 )=disturb_mag;
    disturb_step(disturb_bus(i)+2*num_gen,t_sim<1)=0;
end

fx_pre=@(x) f_x(x,bij_pre,gij_pre,gii_pre,Pm_pre);
Jx_pre=@(x) J_fx(x,bij_pre,gij_pre,gii_pre,Pm_pre);
H=diag([zeros(2*num_gen,1); D_L.^-1; zeros(num_gen,1)]);
[x_sim_cont,t_sim_cont]=tds(x_pre, 0, t_end, del_t, 'trapz', fx_pre, Jx_pre,[],H*disturb_cont);
[x_sim_step,t_sim]=tds(x_pre, 0, t_end, del_t, 'trapz', fx_pre, Jx_pre,[],H*disturb_step);

%%
figure;
subplot(2,2,1); box on; hold on;
plot(t_sim,disturb_cont(disturb_bus+2*num_gen,:),'LineWidth',1.5)
plot(t_sim,disturb_mag*ones(size(t_sim)),'r--','LineWidth',1.5)
plot(t_sim,-disturb_mag*ones(size(t_sim)),'r--','LineWidth',1.5)
y_label=ylabel('$u_1$ ($\Delta P$ in p.u.)'); set(y_label, 'Interpreter', 'latex','FontSize',15,'FontName','Times New Roman');
title_label=title('\textbf{(a) Continuous disturbance}');  set(title_label, 'Interpreter', 'latex','FontSize',15,'FontName','Times New Roman');
xlim([0 t_end]); ylim([-1.5 1.5])
set(gca,'fontsize',15,'FontName','Times New Roman')
subplot(2,2,3); box on; hold on;
plot(t_sim,x_sim_cont(idx_omega,:)/(2*pi),'LineWidth',3.5) % divide by 2*pi to get frequency deviation in Hz
plot(t_sim,0.5*ones(size(t_sim)),'r--','LineWidth',1.5)
plot(t_sim,-0.5*ones(size(t_sim)),'r--','LineWidth',1.5)
y_label=ylabel('$$y_1$ ($\dot{\delta}$ in Hz)'); set(y_label, 'Interpreter', 'latex','FontSize',15,'FontName','Times New Roman');
xlim([0 t_end]); ylim([-1 1])
x_label=xlabel('Time (seconds)'); set(x_label, 'Interpreter', 'latex','FontSize',15,'FontName','Times New Roman');
title_label=title('\textbf{(c) Continuous response}');  set(title_label, 'Interpreter', 'latex','FontSize',15,'FontName','Times New Roman');
set(gca,'fontsize',15,'FontName','Times New Roman')
subplot(2,2,2); box on; hold on;
plot(t_sim,disturb_step(disturb_bus+2*num_gen,:),'Color',[0 0.4470 0.7410],'LineWidth',1.5)
plot(t_sim,disturb_mag*ones(size(t_sim)),'r--','LineWidth',1.5)
plot(t_sim,-disturb_mag*ones(size(t_sim)),'r--','LineWidth',1.5)
y_label=ylabel('$u_2$ ($\Delta P$ in p.u.)'); set(y_label, 'Interpreter', 'latex','FontSize',15,'FontName','Times New Roman');
title_label=title('\textbf{(b) Step disturbance}');  set(title_label, 'Interpreter', 'latex','FontSize',15,'FontName','Times New Roman');
xlim([0 t_end]); ylim([-1.5 1.5])
set(gca,'fontsize',15,'FontName','Times New Roman')
subplot(2,2,4); box on; hold on;
plot(t_sim,x_sim_step(idx_omega,:)/(2*pi),'-','LineWidth',2.5) % divide by 2*pi to get frequency deviation in Hz
plot(t_sim,0.5*ones(size(t_sim)),'r--','LineWidth',1.5)
plot(t_sim,-0.5*ones(size(t_sim)),'r--','LineWidth',1.5)
y_label=ylabel('$y_2$ ($\dot{\delta}$ in Hz)'); set(y_label, 'Interpreter', 'latex','FontSize',15,'FontName','Times New Roman');
xlim([0 t_end]); ylim([-1 1])
x_label=xlabel('Time (seconds)'); set(x_label, 'Interpreter', 'latex','FontSize',15,'FontName','Times New Roman');
title_label=title('\textbf{(d) Step response}');  set(title_label, 'Interpreter', 'latex','FontSize',15,'FontName','Times New Roman');
set(gca,'fontsize',15,'FontName','Times New Roman')

%% Individual plot for presentation
close all;
set(figure,'Position', [50 50 300 250]); hold all; box on;
plot(t_sim,disturb_cont(disturb_bus+2*num_gen,:),'LineWidth',1.5)
plot(t_sim,disturb_mag*ones(size(t_sim)),'r--','LineWidth',1.5)
plot(t_sim,-disturb_mag*ones(size(t_sim)),'r--','LineWidth',1.5)
y_label=ylabel('$u$'); set(y_label, 'Interpreter', 'latex','FontSize',15,'FontName','Times New Roman');
x_label=xlabel('Time'); set(x_label, 'Interpreter', 'latex','FontSize',15,'FontName','Times New Roman');
xlim([0 t_end]); ylim([-1.5 1.5])
set(gca,'fontsize',15,'FontName','Times New Roman')


set(figure,'Position', [50 50 300 250]); hold all; box on;
plot(t_sim,x_sim_step(idx_omega,:)/(2*pi),'LineWidth',3.5) % divide by 2*pi to get frequency deviation in Hz
plot(t_sim,0.5*ones(size(t_sim)),'r--','LineWidth',1.5)
plot(t_sim,-0.5*ones(size(t_sim)),'r--','LineWidth',1.5)
y_label=ylabel('$y$'); set(y_label, 'Interpreter', 'latex','FontSize',15,'FontName','Times New Roman');
x_label=xlabel('Time'); set(x_label, 'Interpreter', 'latex','FontSize',15,'FontName','Times New Roman');
xlim([0 t_end]); ylim([-1 1])
set(gca,'fontsize',15,'FontName','Times New Roman')

%% Angle difference plot
close all
set(figure,'Position', [50 50 300 250]); hold all; box on;
plot(t_sim,disturb_cont(disturb_bus(1)+2*num_gen,:),'r','LineWidth',2)
set(figure,'Position', [50 50 300 250]); hold all; box on;
plot(t_sim,E_N(:,1)'*x_sim_cont(idx_delta,:),'LineWidth',2)
axis([0 30 -0.015 0.02]);
set(figure,'Position', [50 50 300 250]); hold all; box on;
plot(t_sim,E_N(:,11)'*x_sim_cont(idx_delta,:),'LineWidth',2)
%axis([0 30 0.036 0.0385])

%%
close all;
signal=zeros(3,size(t_sim,2));
for i=1:3 
    rand_sig=sum(cos(2*rand(3,1)*t_sim+2*pi*rand(3,1))); 
    signal(i,:)=((1+0.5*(3-i))*rand_sig)/max(abs(rand_sig));
end

set(figure,'Position', [50 50 300 250]); hold all; box on;
plot(t_sim,signal,'LineWidth',1.5)



plot(t_sim,max(abs(signal(1,:))')'*ones(size(t_sim)),'--','Color',[0 0.4470 0.7410],'LineWidth',1.5)
plot(t_sim,max(abs(signal(2,:))')'*ones(size(t_sim)),'--','Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
plot(t_sim,max(abs(signal(3,:))')'*ones(size(t_sim)),'--','Color',[0.9290 0.6940 0.1250],'LineWidth',1.5)
plot(t_sim,-max(abs(signal(1,:))')'*ones(size(t_sim)),'--','Color',[0 0.4470 0.7410],'LineWidth',1.5)
plot(t_sim,-max(abs(signal(2,:))')'*ones(size(t_sim)),'--','Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
plot(t_sim,-max(abs(signal(3,:))')'*ones(size(t_sim)),'--','Color',[0.9290 0.6940 0.1250],'LineWidth',1.5)
%plot(t_sim,min(signal')'*ones(size(t_sim)),'r--')
y_label=ylabel('$u$'); set(y_label, 'Interpreter', 'latex','FontSize',15,'FontName','Times New Roman');
x_label=xlabel('Time'); set(x_label, 'Interpreter', 'latex','FontSize',15,'FontName','Times New Roman');
xlim([0 t_end]); ylim([-3 3])
set(gca,'fontsize',15,'FontName','Times New Roman')
