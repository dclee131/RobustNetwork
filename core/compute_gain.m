clear; close all;
sys_case=39;
compute_gain_method=1; % 1: Direct impulse computation / 2: using / 3: both
loading_level=100;%sys_case=9; compute_gain_method=0;

%% Build 2nd order model with Kron reduction
model_order='2gov'; activate_Lossless=1;
run(['dyn' int2str(sys_case)]) % Get info from data file
run(['model_' model_order]) % Build model

%% Formulate the model
num_nonlin=num_line+num_gen; w_0=E_N'*x_eq(idx_delta);
A=J_fx(x_pre,bij_pre,gij_pre,gii_pre,Pm_pre);
B_u=zeros(num_var,num_gen+num_bus); B_u(idx_pm,1:num_gen)=diag(T_g.^-1); B_u(idx_delta_grid,num_gen+1:end)=diag(D_L.^-1);
B_v=zeros(num_var,num_nonlin); B_v([idx_omega idx_delta_grid],:)=[-diag(M_G.^-1)*E_G*diag(bij_pre); -diag(D_L.^-1)*E_L*diag(bij_pre)];
C_y=zeros(num_gen,num_var); C_y(:,idx_omega)=eye(num_gen);
C_w=zeros(num_nonlin,num_var); C_w(:,idx_delta)=[E_G' E_L'];
%B_u=[zeros(num_gen,num_gen+num_bus); diag(M_G.^-1) zeros(num_gen, num_bus); zeros(num_bus, num_gen) diag(D_L.^-1)];
%B_v=[zeros(num_gen,num_nonlin); -diag(M_G.^-1)*E_G*diag(bij_pre); -diag(D_L.^-1)*E_L*diag(bij_pre)];
%C_y=[zeros(num_gen) eye(num_gen) zeros(num_gen,num_bus)];
%C_w=[E_G' zeros(num_nonlin,num_gen) E_L'];

G{1}=ss(A,B_u,C_y,zeros(num_gen,num_bus+num_gen));
G{2}=ss(A,B_v,C_y,zeros(num_gen,num_nonlin));
G{3}=ss(A,B_u,C_w,zeros(num_nonlin,num_bus+num_gen));
G{4}=ss(A,B_v,C_w,zeros(num_nonlin));

%% Compute Infinity gains
if compute_gain_method==0
    cd(fileparts(which(mfilename)));
    save_file=['../save/gain_mtx_' num2str(sys_case) '_' num2str(loading_level) '.mat'];
    if exist(save_file, 'file') == 2
        load(save_file);
        disp(['Computation time with simulation: ' num2str(cp_sim)])
        disp(['Computation time with reduced order model: ' num2str(cp_reduced)])
    else
        disp('No saved gain data!')
        return;
    end
end

if compute_gain_method==1 || compute_gain_method==3 % compute with direct integration with simulation
    tic;
    for idx_io=1:4
        sys=prescale(G{idx_io});
        for i=1:size(sys,1)
            for j=1:size(sys,2)
                [y_impulse,t_impulse]=impulse(sys(i,j));
                %if abs(y_impulse(end))>1; disp(['Gain not converged!' num2str(y_impulse(end))]); return; end;
                gain_mtx_sim{idx_io}(i,j)=trapz(t_impulse,abs(y_impulse));
            end
        end
    end
    cp_sim=toc;
    disp(['Computation time with simulation: ' num2str(cp_sim)])
end

if compute_gain_method==2 || compute_gain_method==3 % compute with model order reduction
    threshold_reduction=1e-2;
    tic;
    for idx_io=1:4
        [num,den]=tfdata(prescale(G{idx_io}));
        for i=1:size(num,1)
            for j=1:size(num,2)
                [r,p,k] = residue(num{i,j},den{i,j});
                idx_nontriv=find((abs(real(r)./real(p))>threshold_reduction).*(abs(real(p))>1e-5));
                idx_triv=find((abs(real(r)./real(p))<threshold_reduction).*(abs(real(p))>1e-5));
                
                if size(idx_nontriv,1)==0
                    gain_reduced=0;
                else
                    [b_reduced,a_reduced] = residue(r(idx_nontriv),p(idx_nontriv),k);
                    [y_reduced,t_reduced]=impulse(tf(real(b_reduced),real(a_reduced)));
                    gain_reduced=trapz(t_reduced,abs(y_reduced));
                end
                gain_mtx_reduced{idx_io}(i,j)=gain_reduced+sum(abs(real(r(idx_triv))./real(p(idx_triv))));
            end
        end
    end
    cp_reduced=toc;
    disp(['Computation time with reduced order model: ' num2str(cp_reduced)])
end

if compute_gain_method==3
    save_file=['../save/gain_mtx_' num2str(sys_case) '_' num2str(loading_level) '.mat'];
    cd(fileparts(which(mfilename))); save(save_file,'gain_mtx_sim','gain_mtx_reduced','cp_sim','cp_reduced');
    error_min=min(min(gain_mtx_reduced{4}-gain_mtx_sim{4}));
    error_max=max(max(gain_mtx_reduced{4}-gain_mtx_sim{4}));
    disp(['Error range: ' num2str(error_min) '/' num2str(error_max)]);
end

%% Check maximum gain
if 0
    [max_gain,i_max]=max(gain_mtx_sim{4});
    [max_gain,j_max]=max(max_gain);
    i_max=i_max(j_max);
    %i_max=1; j_max=1;
    [y_impulse,t_impulse]=impulse(prescale(G{4}(i_max,j_max)));
    plot(t_impulse,y_impulse)
    disp(['Maximum Gain: ' num2str(max_gain) ' (' num2str(trapz(t_impulse,abs(y_impulse))) ')'])
end