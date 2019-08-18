using MATLAB
using Ipopt
using JuMP

# Setting
sys_case=9;
freq_limit=2*pi*0.5; # 0: no frequency limit, ~0: frequency limit
compute_gain_method=0; # 0: load precomputed gain, 1: direct simulation, 2: model reduction

###### Begin Code #######
# Set up model and compute gain
mat"""
close all;
addpath(genpath(pwd))
sys_case=$sys_case;
compute_gain_method=$compute_gain_method;
compute_gain;
$gain_mtx=gain_mtx_sim;
$w_0=w_0;
$cp_gain=cp_sim;
$line_idx=line_frto+num_gen;
"""
size_w=size(w_0,1);
size_u=size(gain_mtx[1],2);
size_y=size(gain_mtx[1],1);

run_c=1:size_u;
#run_c=(size_y+1:size_u)';
#run_c=round.(Int8, line_idx);
#run_c=[1 9 16; 3 15 27]+size_y;

u_bar=zeros(size_u); y_bar=zeros(size_y); w_bar=zeros(size_w); sinw_bar=zeros(size_w);
u_bar_save=zeros(size(run_c,1));
for i=1:size(run_c,1)
 #m = Model(solver=IpoptSolver(print_level=0));
 m = Model(solver=IpoptSolver());
 @variable(m, 0 <= w[1:size_w] <= pi/2, start=0)
 @variable(m, 0 <= sinw[1:size_w] <= 1, start=0)
 @variable(m, 0 <= u[1:size_u] <= 100, start=0)
 @variable(m, mu, start=0)

 @NLconstraint(m, myconstr[i=1:size_w], sinw[i] <= sin(abs.(w_0[i])+w[i]))
 @constraint(m,gain_mtx[3]*u.<=w-gain_mtx[4]*(diagm(cos.(w_0))*w-sinw+sin.(abs.(w_0))))
 if freq_limit>0
  @constraint(m,gain_mtx[1]*u+gain_mtx[2]*(diagm(cos.(w_0))*w-sinw+sin.(abs.(w_0))).<=freq_limit)
 end
 #c=zeros(size_u); c[run_c[i,:]]=1;
 @constraint(m,mu*ones(size(run_c,2)).<=u[run_c[i,:]])
 @objective(m, Max, mu)
 status=solve(m);
 u_bar=getvalue(u);
 w_bar=getvalue(w);

 sinw_bar=getvalue(sinw);
 y_bar=gain_mtx[1]*u_bar+gain_mtx[2]*(diagm(cos.(w_0))*w_bar-sinw_bar+sin.(abs.(w_0)))
 u_bar_save[i,:]=getvalue(mu);
end

println("w = ", w_bar); println("y = ", y_bar);
println("Result Summary")
println("Maximum perturbation: ", u_bar_save);
println("Gain Computation Time: ", cp_gain, " seconds")
#println("IPOPT Optimization Time: ", cp_gain, " seconds")

# Visualize result on the network
if sum(sys_case.==[9 14 39])>0
 if size(run_c,2)==1
  mat"VisualizeNetwork(sys_case,$u_bar_save,$w_bar,[0.5 1],[0 10],0);"
  #mat" VisualizeNetwork(sys_case,$u_bar_save,$w_bar,[0.2 5],[0 100],0);"
 end
end

sgt_mtx=eye(size_w)-gain_mtx[4]*diagm(cos.(w_0)+(sin.(abs.(w_0))-sinw_bar)./w_bar);
#mat"figure; imagesc(gain_mtx{4}); colorbar;" # Gain matrix visualization
#mat"figure; imagesc(sgt_mtx); colorbar;" # M matrix visualization
#mat"figure; scatter(bij_pre,$w_bar)" # Relationship between line impedence and angle bound
#mat"figure; scatter(diag(E_N*diag(bij_pre)*E_N'),$u_bar_save)" # Bus transfer capacity vs perturbation bound
