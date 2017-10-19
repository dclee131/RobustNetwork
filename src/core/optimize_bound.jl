using MATLAB
using Ipopt
using JuMP

# Setting
sys_case=9;
pert_bus=6;

###### Begin Code #######
# Set up model and compute gain
mat"""
addpath(genpath(pwd))
sys_case=$sys_case;
method=1;
compute_gain;
$gain_mtx=gain_mtx;
$w_0=E*delta0;
$num_line=num_line;
$num_bus=num_bus;
$idx_gen=idx_gen;
"""
size_w=convert(Int64,num_line);
size_u=convert(Int64,num_bus);

# Formulate and solve optimization problem
c=zeros(size_u); c[pert_bus]=1;

m = Model(solver=IpoptSolver());
@variable(m, 0 <= w[1:size_w] <= pi/2,start=0)
@variable(m, 0 <= sinw[1:size_w] <= 1,start=0)
@variable(m, 0 <= u[1:size_u] <= pi/2,start=0)

@NLconstraint(m, myconstr[i=1:size_w], sinw[i] <= sin(abs.(w_0[i])+w[i]))
@constraint(m,gain_mtx[3]*u.<=(eye(size_w)-gain_mtx[4]*diagm(cos.(w_0)))*w-gain_mtx[4]*sin.(abs.(w_0))+gain_mtx[4]*sinw)
@objective(m, Max, c'*u)
solve(m)
w_bar=getvalue(w); u_bar=getvalue(u);
println("w = ", w_bar); println("u = ", u_bar);

# Visualize result on the network
mat"""
if sum(sys_case==[9 39])
 close all; VisualizeNetwork(sys_case,$u_bar,$w_bar);
end
"""
