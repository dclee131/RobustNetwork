clear; close all;
sys_case=9; loading_level=100;
step_mag=0.5; pert_bus=2;

% sys_case=9; loading_level=100; step_mag=0.5; pert_bus=2; % Frequency plot setting

%% Build 2nd order model with Kron reduction
model_order='2gov'; activate_Lossless=1;
run(['dyn' int2str(sys_case)]) % Get info from data file
run(['model_' model_order]) % Build model
step_change=zeros(num_bus,1); step_change(pert_bus)=step_mag;

%% Run simulation
method='trapz';
del_t=0.05; % times step for tds
t_fault=1; % ticme that the fault is applied
t_end=10; % end time for time domain simulation (tds)

fx_pre=@(x) f_x(x,bij_pre,gij_pre,gii_pre,Pm_pre);
Jx_pre=@(x) J_fx(x,bij_pre,gij_pre,gii_pre,Pm_pre);
[x_sim_pre,t_sim_pre]=tds(x_pre, 0, t_fault, del_t, method, fx_pre, Jx_pre); % pre-contingency simulation

fx_post=@(x) f_x(x,bij_post,gij_post,gii_post+step_change,Pm_pre);
Jx_post=@(x) J_fx(x,bij_post,gij_post,gii_post+step_change,Pm_pre);
[x_sim_post,t_sim_post]=tds(x_sim_pre(:,end), t_fault, t_end, del_t, method, fx_post, Jx_post); % post-contingency simulation

x_sim=[x_sim_pre x_sim_post];
t_sim=[t_sim_pre t_sim_post];

%% Plot simulation
figure;
subplot(2,1,1); hold on; grid on; box on; plot(t_sim,x_sim(idx_delta,:))
set(gca,'FontSize',15,'FontName','Times New Roman'); ylabel('\delta (rad)'); xlim([0 t_end]);
subplot(2,1,2); hold on; grid on; box on; plot(t_sim,x_sim(idx_omega,:)/2/pi+60) % divide by 2*pi to get frequency deviation in Hz
set(gca,'FontSize',15,'FontName','Times New Roman'); ylabel('\omega (Hz)'); xlim([0 t_end]);

