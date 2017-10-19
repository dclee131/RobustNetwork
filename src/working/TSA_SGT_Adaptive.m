clear; close all;
sopt_gurobi=sdpsettings('verbose',0,'solver','gurobi');

%% Settings
sys_case=9; % IEEE system

run(['dyn' int2str(sys_case)]) % Get info from data file
slack_bus=SW.con(1);
num_bus=size(Bus.con,1);
num_line=size(Line.con,1);
num_gen=size(Syn.con,1);
num_load=num_bus-num_gen;

eye_bus=eye(num_bus); idx_gen=setdiff(Syn.con(:,1),slack_bus); idx_load=setdiff(Bus.con(:,1),Syn.con(:,1));

Syn.con(:,19)=2*ones(size(Syn.con(:,19))); % Add damping
Syn.con(:,4)=10;
M=Syn.con(:,18)./(2*pi*Syn.con(:,4));
D_G=Syn.con(:,19)./(2*pi*Syn.con(:,4));
D_L=30*ones(num_load,1)./(2*pi*Syn.con(1,4));
v_gen=[SW.con(:,4); PV.con(:,5)];
line_frto=Line.con(:,1:2);
Z=Line.con(:,8)+1i*Line.con(:,9);
E=eye_bus(line_frto(:,1),:)-eye_bus(line_frto(:,2),:);
Y=E'*diag(Z.^-1)*E;
B_line=-imag(1./Z);

Pgen=zeros(num_bus,1); Pgen(PV.con(:,1))=PV.con(:,4);
Sload=zeros(num_bus,1); Sload(PQ.con(:,1))=(PQ.con(:,4)+1i*PQ.con(:,5));

% Solve for Equilibrium
x_eq=NR_ss(Y,Pgen-Sload,idx_load,v_gen,slack_bus);
delta_eq=x_eq(1:num_bus);
V_eq=x_eq(num_bus+1:end).*(cos(x_eq(1:num_bus))+1i*sin(x_eq(1:num_bus)));
I_eq=Y*V_eq;
S_inj=V_eq.*conj(I_eq);
P_inj=real(S_inj);
P_G=P_inj(idx_gen);
P_L=P_inj(idx_load);
%P_inj=E'*diag(B_line)*sin(E*delta_eq);

num_gen=num_gen-1; num_bus=num_bus-1; idx_bus=setdiff(Bus.con(:,1),slack_bus);
E_G=E(:,idx_gen); E_L=E(:,idx_load); E=[E_G E_L];
M=M(Syn.con(:,1)~=slack_bus); D_G=D_G(Syn.con(:,1)~=slack_bus);
idx_delta_gen=1:num_gen; idx_omega=num_gen+1:2*num_gen; idx_delta_load=2*num_gen+1:num_bus+num_gen;
idx_delta=[1:num_gen 2*num_gen+1:num_bus+num_gen];

%% Set up simulation model
f_sys= @(x) [x(idx_omega);
    diag(1./M)*(-D_G.*x(idx_omega)-E_G'*diag(B_line)*sin(E*x(idx_delta))+P_G);
    diag(1./D_L)*(-E_L'*diag(B_line)*sin(E*x(idx_delta))+P_L)];
J_sys=@(x) [zeros(num_gen) eye(num_gen) zeros(num_gen,num_load);
    -diag(1./M)*E_G'*diag(B_line)*diag(cos(E*x(idx_delta)))*E_G -diag(D_G./M) -diag(1./M)*E_G'*diag(B_line)*diag(cos(E*x(idx_delta)))*E_L;
    -diag(1./D_L)*E_L'*diag(B_line)*diag(cos(E*x(idx_delta)))*E_G zeros(num_load,num_gen) -diag(1./D_L)*E_L'*diag(B_line)*diag(cos(E*x(idx_delta)))*E_L];

delta0=(E'*diag(B_line)*E)\[P_G; P_L];
x0=[delta0(1:num_gen); zeros(num_gen,1); delta0(num_gen+1:num_bus)];

% Jacobian check
x=x0; J_check=zeros(size(x));
for i=1:size(x,1)
    x_unit=zeros(size(x,1),1); x_unit(i)=1;
    J_check(:,i)=(f_sys(x+0.001*x_unit)-f_sys(x))/0.001;
end
J_fail=max(max(abs(J_sys(x)-J_check)));

x_eq = NR(f_sys,J_sys,x0);
if max(x_eq)>1e3; return; end

%% Formulate the model
A=[zeros(num_gen) eye(num_gen) zeros(num_gen,num_load);
    -diag(M.^-1)*E_G'*diag(B_line)*diag(cos(E*x_eq(idx_delta)))*E_G -diag(D_G./M) -diag(M.^-1)*E_G'*diag(B_line)*diag(cos(E*x_eq(idx_delta)))*E_L;
    -diag(D_L.^-1)*E_L'*diag(B_line)*diag(cos(E*x_eq(idx_delta)))*E_G zeros(num_load,num_gen) -diag(D_L.^-1)*E_L'*diag(B_line)*diag(cos(E*x_eq(idx_delta)))*E_L];
B=[zeros(num_gen,num_bus); diag(M.^-1) zeros(num_gen, num_load); zeros(num_load, num_gen) diag(D_L.^-1)];
B_v=[zeros(num_gen,num_line); -diag(M.^-1)*E_G'*diag(B_line); -diag(D_L.^-1)*E_L'*diag(B_line)];

G{1}=ss(A,B,eye(num_bus+num_gen),zeros(num_bus+num_gen,num_bus));
G{2}=ss(A,B_v,eye(num_bus+num_gen),zeros(num_bus+num_gen,num_line));
G{3}=ss(A,B,[E_G zeros(num_line,num_gen) E_L],zeros(num_line,num_bus));
G{4}=ss(A,B_v,[E_G zeros(num_line,num_gen) E_L],zeros(num_line));

%% Compute Infinity gains
for idx_io=1:4
    sys=G{idx_io};
    for i=1:size(sys,1)
        for j=1:size(sys,2)
            [y_impulse,t_impulse]=impulse(sys(i,j));
            gain_mtx{idx_io}(i,j)=trapz(t_impulse,abs(y_impulse));
        end
    end
    disp(['Progress: ' num2str(idx_io) '/4'])
end




%%
%bounding_vector=ones(num_bus,1); % Uniform
%bounding_vector=zeros(num_bus,1); bounding_vector(2)=1; % Select a bus
%bounding_vector=zeros(num_bus+1,1); bounding_vector(idx_gen)=1; bounding_vector(slack_bus)=[]; % Select Generators
bounding_vector=zeros(num_bus+1,1); bounding_vector(idx_load)=1; bounding_vector(slack_bus)=[]; % Select Loads

w_max_plot=0.001:0.01:0.5; u_max=[];
for i=1:size(w_max_plot,2) % choose bound on delta
    w_max=w_max_plot(i);
    uncertainty_up=abs(cos(E*x_eq(idx_delta))-(sin(E*x_eq(idx_delta)+w_max)-sin(E*x_eq(idx_delta)))/w_max);
    uncertainty_down=abs(cos(E*x_eq(idx_delta))-(sin(E*x_eq(idx_delta))-sin(E*x_eq(idx_delta)-w_max))/w_max);
    gain_nonlinearity=diag(max(uncertainty_up,uncertainty_down));
    u_max(i)=min((gain_mtx{3}*bounding_vector).^-1.*(eye(num_line)-gain_mtx{4}*gain_nonlinearity)*w_max*ones(num_line,1));
end

close all;
figure; hold all; box on; grid on;
plot(w_max_plot,u_max)

%% Direct Simulation
if 0
    del_t_sim=0.1; t_end=200;
    t_sim=0:del_t_sim:t_end;
    P_rng=0.1:0.2:5;
    u_intensity=bounding_vector*P_rng; %u_intensity=u_intensity*0.1;
    x_sim=zeros(size(x_eq,1),t_end/del_t_sim);
    delta_max_sim=inf*ones(size(P_rng));
    
    for i=1:size(u_intensity,2)
        % Simulate worst case
        x_sim(:,1)=x_eq;
        for t_idx=1:t_end/del_t_sim
            f=@(x) -x+x_sim(:,t_idx)+del_t_sim/2*(f_sys(x_sim(:,t_idx))+B*u_intensity(:,i)+f_sys(x)+B*u_intensity(:,i));
            J=@(x) -eye(num_gen+num_bus)+del_t_sim/2*J_sys(x);
            x_sim(:,t_idx+1)=NR(f,J,x_sim(:,t_idx));
            if max(x_sim(:,t_idx))>1e3; break; end
        end
        delta_max_sim(i)=max(max(abs(E*(x_sim(idx_delta,:)-x_eq(idx_delta)))));
        if delta_max_sim(i)>10; break; end;
        %figure; plot(t_sim,x_sim);
    end
    plot([0 delta_max_sim+1e10*(delta_max_sim>pi)],[0 P_rng]); %xlim([0 pi])
    set(gca,'FontSize',15,'FontName','Times New Roman'); xlabel('||w||_\infty'); ylabel('\Delta P_{max}');
    legend('Lower bound','Upper bound','Orientation','horizontal')
    xlim([0 pi])
end

%% Frequency deviation bound
[u_max_max,idx_u_max]=max(u_max);
gain_mtx{1}*u_max_max*bounding_vector+gain_mtx{2}*gain_nonlinearity*w_max_plot(idx_u_max)*ones(num_line,1);
w_max_plot(idx_u_max)

%%
c=bounding_vector;
w_max=pi/4*ones(num_line,1); %initialize

idx_w=1;
w_opt=w_max;
w_i_min=0.0001; w_i_max=pi/3;
for i=1:10
    w_opt(idx_w)=(w_i_min+w_i_max)/2;
    w_opt(idx_w)=0.1*i;
    %w_opt=w_max_plot(idx_u_max)*ones(num_line,1);
    uncertainty_up=cos(E*x_eq(idx_delta)).*w_opt-(sin(E*x_eq(idx_delta)+w_opt)-sin(E*x_eq(idx_delta)));
    uncertainty_down=cos(E*x_eq(idx_delta)).*w_opt-(sin(E*x_eq(idx_delta))-sin(E*x_eq(idx_delta)-w_opt));
    gain_nonlinearity=max(abs(uncertainty_up),abs(uncertainty_down));
    u=sdpvar(num_bus,1);
    constr=[gain_mtx{3}*u<=w_opt-gain_mtx{4}*gain_nonlinearity; u>=0];
    result=optimize(constr,-c'*u,sopt_gurobi);
    
    [idle,idx_binding]=min(w_opt-gain_mtx{4}*gain_nonlinearity-gain_mtx{3}*value(u));
    % Find gradient
    
    gradient_nonlinearity=cos(E*x_eq(idx_delta)+w_opt)-cos(E*x_eq(idx_delta));
    gradient_nonlinearity=gradient_nonlinearity.*sign(uncertainty_up);
    full_grad(abs(uncertainty_up)<abs(uncertainty_down))=gradient_nonlinearity(abs(uncertainty_up)<abs(uncertainty_down));
    
    gradient_nonlinearity=cos(E*x_eq(idx_delta))-cos(E*x_eq(idx_delta)-w_opt);
    gradient_nonlinearity=gradient_nonlinearity.*sign(uncertainty_down);
    full_grad(abs(uncertainty_up)>abs(uncertainty_down))=gradient_nonlinearity(abs(uncertainty_up)>abs(uncertainty_down));
    
    gradient=1-gain_mtx{4}(idx_binding,idx_w)*full_grad(idx_w);
    %if gradient<0 || min(min(inv(eye(num_line)-gain_mtx{4}*diag(gain_nonlinearity./w_opt))))<0; 
    if gradient<0
        w_i_max=w_opt(idx_w); 
    else
        w_i_min=w_opt(idx_w);
    end
    d(i)=value(c'*u);
    dd(i)=min(w_opt-gain_mtx{4}*gain_nonlinearity);
    ddd(i)=gradient;
    z(i)=gain_nonlinearity(idx_w);
    zz(i)=full_grad(idx_w);
    q(i)=idx_binding;
end
value(c'*u)
w=(w_i_max+w_i_min)/2;
disp([num2str(pi/3) ' ' num2str(w)])
close all; subplot(2,1,1); hold all; plot(dd,'r'); plot(ddd,'b--');
subplot(2,1,2); hold all; plot(z,'r'); plot(zz,'b--');