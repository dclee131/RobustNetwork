clear; close all;

%% Settings
sys_case=2; % IEEE system

%% Set up variables
run(['dyn' int2str(sys_case)]) % Get info from data file
slack_bus=SW.con(1);
num_bus=size(Bus.con,1);
num_line=size(Line.con,1);
num_gen=size(Syn.con,1);
num_load=num_bus-num_gen;

eye_bus=eye(num_bus); idx_gen=Syn.con(:,1); idx_load=setdiff(1:num_bus,idx_gen);
idx_delta=1:num_gen; idx_omega=num_gen+1:2*num_gen; num_var=2*num_gen;

M_gen=Syn.con(:,18)./(2*pi*Syn.con(:,4));
M_T=sum(M_gen);
Syn.con(:,19)=2*ones(size(Syn.con(:,19))); % Add damping
D_gen=Syn.con(:,19)./(2*pi*Syn.con(:,4));
v_gen=[SW.con(:,4); PV.con(:,5)];
xd_p=Syn.con(:,9);
line_frto=Line.con(:,1:2);
Z_line=Line.con(:,8)+1i*Line.con(:,9);
E=eye_bus(line_frto(:,1),:)-eye_bus(line_frto(:,2),:);
Y=E'*diag(Z_line.^-1)*E;

Pgen=zeros(num_bus,1); Pgen(PV.con(:,1))=PV.con(:,4);
Sload=zeros(num_bus,1); Sload(PQ.con(:,1))=(PQ.con(:,4)+1i*PQ.con(:,5));

% Convert to constant impedence load
x_eq=NR_ss(Y,Pgen-Sload,idx_load,v_gen,slack_bus);
V_eq=x_eq(num_bus+1:end).*(cos(x_eq(1:num_bus))+1i*sin(x_eq(1:num_bus)));
I_eq=Y*V_eq;
S_inj=V_eq.*conj(I_eq);
y_load=conj(Sload)./V_eq.^2;
YN_pre=E'*diag(Z_line.^-1)*E+diag(y_load);

% Include stator impedence for network reduction
Y_pre=zeros(num_bus+num_gen);
Y_pre([1:num_gen,num_gen+idx_gen'],[1:num_gen,num_gen+idx_gen'])=[diag((1i*xd_p).^-1) diag(-(1i*xd_p).^-1); diag(-(1i*xd_p).^-1) diag((1i*xd_p).^-1)];
Y_pre(num_gen+1:end,num_gen+1:end)=Y_pre(num_gen+1:end,num_gen+1:end)+YN_pre;

% Pre-contingency Equilibrium
x_eq_pre=NR_ss(YN_pre,Pgen,idx_load,v_gen,slack_bus);
V_eq_pre=x_eq_pre(num_bus+1:end).*(cos(x_eq_pre(1:num_bus))+1i*sin(x_eq_pre(1:num_bus)));
I_eq_pre=YN_pre*V_eq_pre;
Pgen_pre=real(V_eq_pre(idx_gen).*conj(I_eq_pre(idx_gen)));
Eeq_pre=abs(V_eq_pre(idx_gen)+1i*xd_p.*I_eq_pre(idx_gen));
delta_eq_pre=angle(V_eq_pre(idx_gen)+1i*xd_p.*I_eq_pre(idx_gen));
x_eq_pre=[delta_eq_pre-M_gen'*delta_eq_pre/M_T; zeros(num_gen,1)];

% Apply Kron Reduction
Y_pre_kron=Y_pre(1:num_gen,1:num_gen)-Y_pre(1:num_gen,num_gen+1:end)*(Y_pre(num_gen+1:end,num_gen+1:end)\Y_pre(num_gen+1:end,1:num_gen));
edge_kron=nchoosek(1:num_gen,2);
E_kron=zeros(size(edge_kron,1),num_gen);
for i=1:size(E_kron,1)
    E_kron(i,edge_kron(i,1))=1;  E_kron(i,edge_kron(i,2))=-1;
    gij_pre_kron(i,1)=real(Y_pre_kron(E_kron(i,:)==1,E_kron(i,:)==-1))*Eeq_pre(E_kron(i,:)==1)*Eeq_pre(E_kron(i,:)==-1);
    bij_pre_kron(i,1)=imag(Y_pre_kron(E_kron(i,:)==1,E_kron(i,:)==-1))*Eeq_pre(E_kron(i,:)==1)*Eeq_pre(E_kron(i,:)==-1);
end
P_pre_kron=Pgen_pre-real(diag(Y_pre_kron)).*Eeq_pre.^2;

 %% Model
f_pre_kron= @(x) [x(idx_omega)-1/M_T*M_gen'*x(idx_omega); diag(1./M_gen)*(-D_gen.*x(idx_omega)-E_kron'*diag(bij_pre_kron)*sin(E_kron*x(idx_delta))-abs(E_kron)'*diag(gij_pre_kron)*cos(E_kron*x(idx_delta))+P_pre_kron)-1/M_T*(sum(P_pre_kron)-2*ones(size(E_kron'))*diag(gij_pre_kron)*cos(E_kron*x(idx_delta)))];
J_pre_kron=@(x) [zeros(num_gen) eye(num_gen)-ones(num_gen,1)*M_gen'/M_T; diag(1./M_gen)*(-E_kron'*diag(bij_pre_kron)*diag(cos(E_kron*x(idx_delta)))*E_kron+abs(E_kron)'*diag(gij_pre_kron)*diag(sin(E_kron*x(idx_delta)))*E_kron)-2/M_T*ones(size(E_kron'))*diag(gij_pre_kron)*diag(sin(E_kron*x(idx_delta)))*E_kron -diag(D_gen./M_gen)];

%% Linearize and get peak
J_SGT=@(x,w) [zeros(num_gen) eye(num_gen)-ones(num_gen,1)*M_gen'/M_T; diag(1./M_gen)*(-E_kron'*diag(bij_pre_kron)*diag(cos(E_kron*x(idx_delta)))*E_kron+abs(E_kron)'*diag(gij_pre_kron)*diag(sin(E_kron*x(idx_delta)))*E_kron)-2/M_T*ones(size(E_kron'))*diag(gij_pre_kron)*diag(sin(E_kron*x(idx_delta)))*E_kron -diag(D_gen./M_gen)];
A=J_pre_kron(x_eq_pre);
B=[zeros(num_var/2); eye(num_var/2)];
C=[eye(num_var/2) zeros(num_var/2)];
D=zeros(num_var/2);
sys=ss(A,B,C,D);
for i=1:num_var/2
    for j=1:num_var/2
        tf_mtx(i,j)=getPeakGain(sys(i,j));
    end
end

%%
w_max=0.3;% choose bound on delta
uncertainty_up=abs(cos(x_eq_pre(idx_delta))-(sin(x_eq_pre(idx_delta)+w_max)-sin(x_eq_pre(idx_delta)))/w_max);
uncertainty_down=abs(cos(x_eq_pre(idx_delta))-(sin(x_eq_pre(idx_delta))-sin(x_eq_pre(idx_delta)-w_max))/w_max);
max(uncertainty_up,uncertainty_down)


