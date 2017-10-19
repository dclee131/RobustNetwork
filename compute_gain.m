%% Settings
%clear; close all; sys_case=39; method=1;
%method=1; % 1: Direct impulse computation / 2: using / 3: both

%%
if ~exist('sys_case','var'); sys_case=9; end % default IEEE system
if ~exist('method','var'); method=1; end % default IEEE system

run(['dyn' int2str(sys_case)]) % Get info from data file
slack_bus=SW.con(1);
num_bus=size(Bus.con,1);
num_line=size(Line.con,1);
num_gen=size(Syn.con,1);
num_load=num_bus-num_gen;
num_var=num_bus+num_gen;

eye_bus=eye(num_bus); idx_gen=Syn.con(:,1); idx_load=setdiff(Bus.con(:,1),Syn.con(:,1));

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
P_G=P_inj(idx_gen); P_G(Syn.con(:,1)==slack_bus)=-sum(P_inj(setdiff(Bus.con(:,1),slack_bus)));
P_L=P_inj(idx_load);
%P_inj=E'*diag(B_line)*sin(E*delta_eq);

%num_gen=num_gen-1; num_bus=num_bus-1; idx_bus=setdiff(Bus.con(:,1),slack_bus);
%M=M(Syn.con(:,1)~=slack_bus); D_G=D_G(Syn.con(:,1)~=slack_bus);
E_G=E(:,idx_gen); E_L=E(:,idx_load); E=[E_G E_L];
idx_delta_gen=1:num_gen; idx_omega=num_gen+1:2*num_gen; idx_delta_load=2*num_gen+1:num_bus+num_gen;
idx_delta=[1:num_gen 2*num_gen+1:num_bus+num_gen];

%% Set up simulation model
f_sys= @(x) [x(idx_omega);
    diag(1./M)*(-D_G.*x(idx_omega)-E_G'*diag(B_line)*sin(E*x(idx_delta))+P_G);
    diag(1./D_L)*(-E_L'*diag(B_line)*sin(E*x(idx_delta))+P_L)];
J_sys=@(x) [zeros(num_gen) eye(num_gen) zeros(num_gen,num_load);
    -diag(1./M)*E_G'*diag(B_line)*diag(cos(E*x(idx_delta)))*E_G -diag(D_G./M) -diag(1./M)*E_G'*diag(B_line)*diag(cos(E*x(idx_delta)))*E_L;
    -diag(1./D_L)*E_L'*diag(B_line)*diag(cos(E*x(idx_delta)))*E_G zeros(num_load,num_gen) -diag(1./D_L)*E_L'*diag(B_line)*diag(cos(E*x(idx_delta)))*E_L];

delta0=[0; (E(:,2:end)'*diag(B_line)*E(:,2:end))\[P_G(2:end); P_L]];
x0=[delta0(1:num_gen); zeros(num_gen,1); delta0(num_gen+1:num_bus)];

% Jacobian check
x=x0; J_check=zeros(size(x));
for i=1:size(x,1)
    x_unit=zeros(size(x,1),1); x_unit(i)=1;
    J_check(:,i)=(f_sys(x+0.001*x_unit)-f_sys(x))/0.001;
end
J_fail=max(max(abs(J_sys(x)-J_check)));

f_eq=@(x) [f_sys(x); x(1)];
J_eq=@(x) [J_sys(x); 1 zeros(1,num_var-1)];
x_eq = NR(f_eq,J_eq,x0);
if max(x_eq)>1e3; return; end

%% Formulate the model
A=[zeros(num_gen) eye(num_gen) zeros(num_gen,num_load);
    -diag(M.^-1)*E_G'*diag(B_line)*diag(cos(E*x_eq(idx_delta)))*E_G -diag(D_G./M) -diag(M.^-1)*E_G'*diag(B_line)*diag(cos(E*x_eq(idx_delta)))*E_L;
    -diag(D_L.^-1)*E_L'*diag(B_line)*diag(cos(E*x_eq(idx_delta)))*E_G zeros(num_load,num_gen) -diag(D_L.^-1)*E_L'*diag(B_line)*diag(cos(E*x_eq(idx_delta)))*E_L];
B_u=[zeros(num_gen,num_bus); diag(M.^-1) zeros(num_gen, num_load); zeros(num_load, num_gen) diag(D_L.^-1)];
B_v=[zeros(num_gen,num_line); -diag(M.^-1)*E_G'*diag(B_line); -diag(D_L.^-1)*E_L'*diag(B_line)];
C_y=[zeros(num_gen) eye(num_gen) zeros(num_gen,num_load)];
C_w=[E_G zeros(num_line,num_gen) E_L];

G{1}=ss(A,B_u,C_y,zeros(num_gen,num_bus));
G{2}=ss(A,B_v,C_y,zeros(num_gen,num_line));
G{3}=ss(A,B_u,C_w,zeros(num_line,num_bus));
G{4}=ss(A,B_v,C_w,zeros(num_line));

%% Compute Infinity gains
if method==1 || method==3 % compute with direct integration with simulation
    tic;
    for idx_io=1:4
        sys=prescale(G{idx_io});
        for i=1:size(sys,1)
            for j=1:size(sys,2)
                [num,den] = tfdata(sys(i,j));
                [r,p,k] = residue(num{1},den{1});
                [y_impulse,t_impulse]=impulse(sys(i,j));
                %if abs(y_impulse(end))>1; disp(['Gain not converged!' num2str(y_impulse(end))]); return; end;
                gain_mtx_sim{idx_io}(i,j)=trapz(t_impulse,abs(y_impulse));
            end
        end
        %disp(['Progress: ' num2str(idx_io) '/4'])
    end
    gain_mtx=gain_mtx_sim;
    cp_direct=toc;
    disp(['Computation time for direct approach with simulation: ' num2str(cp_direct)])
end

if method==2 || method==3 % compute with direct integration with analytical expression
    tic;
    for idx_io=1:4
        sys=prescale(G{idx_io});
        for i=1:size(sys,1)
            for j=1:size(sys,2)
                [num,den] = tfdata(sys(i,j));
                [r,p,k] = residue(num{1},den{1});
                idx_nontriv=find((abs(real(r)./real(p))>1e-3).*(abs(real(p))>1e-5));
                idx_triv=find((abs(real(r)./real(p))<1e-3).*(abs(real(p))>1e-5));
                t_end=max(abs(log(real(abs(r(idx_nontriv)))/1e-5)./real(p(idx_nontriv))));
                time_step=min(2*pi/10/max(abs(imag(p(idx_nontriv)))),t_end/2);
                t_impulse=0:time_step:t_end;
                y_impulse=real(r(idx_nontriv).'*exp(p(idx_nontriv)*t_impulse));
                %if abs(y_impulse(end))>1; disp(['Gain not converged!' num2str(y_impulse(end))]); return; end;
                gain_mtx_anly{idx_io}(i,j)=trapz(t_impulse,abs(y_impulse))+sum(abs(real(r(idx_triv))./real(p(idx_triv))));
            end
        end
        %disp(['Progress: ' num2str(idx_io) '/4'])
    end
    gain_mtx=gain_mtx_anly;
    cp_direct=toc;
    disp(['Computation time for direct approach with analytical expression: ' num2str(cp_direct)])
end

if method==4
    tic;
    [num,den] = tfdata(G{1}(1,1));
    sys_order=size(den{1},2);
    [r,p,k] = residue(num{1},den{1});
    figure; hold all;
    for idx_order=1:sys_order
        num=zeros(1,sys_order); num(idx_order)=1;
        [unit_impulse{idx_order},t_impulse{idx_order}]=impulse(tf(num,den{1}));
        plot(t_impulse{idx_order},unit_impulse{idx_order})
    end
    
    for idx_io=1:4
        sys=G{idx_io};
        [num,den] = tfdata(sys);
        for i=1:size(sys,1)
            for j=1:size(sys,2)
                y_impulse=sum(diag(num{i,j})*unit_impulse);
                gain_mtx{idx_io}(i,j)=trapz(t_impulse,abs(y_impulse));
            end
        end
    end
    cp_improv=toc;
    disp(['Computation time for direct approach: ' num2str(cp_improv)])
end
if method==3; max(max(abs((gain_mtx_anly{4}-gain_mtx_sim{4})))); end

% Test eigenvalue w.r.t. system size
%n_rng=1:100; for n=n_rng; eig_size(n)=min(real(eig(eye(n)-0.1*ones(n,n)))); end; plot(n_rng,eig_size)

%theta_eq=E*delta_eq;
%diag(cos(theta_eq)+(sin(abs(theta_eq))-sin(abs(theta_eq)+w_bar))/w_bar)
%(eye(size_w)-gain_mtx[4]*diagm(cos.(w_0)))*w-gain_mtx[4]*sin.(abs(w_0))+gain_mtx[4]*sinw)
%eye(num_line)-gain_mtx{4}
%imagesc(gain_mtx{4}); colorbar;
