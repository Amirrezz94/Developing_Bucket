using JuMP
using Ipopt
using Plots
##Supplementary Files
include("CollMat.jl")
##Required Parameters
Nx = 2;
Nu = 1;
NFE = 20;
NCP = 4;
x0= [0.0, 1.0];
dx0=[0,0];
u0=0;
q0=0;
dq0=0;
p = 1;
M = Collocation_Matrix(NCP)
dt = 1.0;
##Model Defining
m1 = Model(Ipopt.Optimizer)
##===========================================================
@variable(m1,  x[1:Nx, 1:NFE, 1:NCP+1])
@variable(m1, dx[1:Nx, 1:NFE, 1:NCP])
@variable(m1,  u[1:Nu, 1:NFE       ])
@variable(m1,  q[1,    1:NFE, 1:NCP])
@variable(m1, dq[1,    1:NFE, 1:NCP])
##===========================================================
for nx in 1:Nx, nfe in 1:NFE, ncp in 1:NCP, nu in 1:Nu
    #Bounds
    set_lower_bound(x[nx, nfe, ncp], -10)
    set_upper_bound(x[nx , nfe, ncp], 10)
    set_lower_bound(u[nu, nfe], -2)
    set_upper_bound(u[nu, nfe], 2)
    #Start Values
    set_start_value(x[nx, nfe, ncp],    x0[nx])
    set_start_value(dx[nx, nfe, ncp], dx0[nx])
    set_start_value(u[nfe],         u0[nu])
    set_start_value(q[1, nfe, ncp],     q0)
    set_start_value(dq[1, nfe, ncp],    dq0)
end
##Set the ODEs===========================================================
@NLconstraints(m1, begin
    ODE1[nfe in 1:NFE, ncp in 1:NCP],    dx[1, nfe, ncp]      == (p - x[2, nfe, ncp] ^2) * x[1, nfe, ncp] - x[2, nfe, ncp] + u[nfe];
    ODE2[nfe in 1:NFE, ncp in 1:NCP],    dx[2, nfe, ncp]      == x[1, nfe, ncp];
    #ODEn[nfe in 1:NFE, ncp in 1:NCP],   ...
end)
##Quadrature===========================================================
@NLconstraints(m1, begin
      Quad_Diff[nfe in 1:NFE, ncp in 1:NCP], dq[1, nfe, ncp]  == (x[1, nfe, ncp])^2 + (x[2, nfe, ncp])^2 + (u[nfe])^2
end)
##Collcation!===========================================================
@NLconstraints(m1, begin
        Coll_Eq_Diff[    nx in 1:Nx,  nfe in 1:NFE,    ncp in 1:NCP],     x[nx, nfe, ncp+1] == x[nx, nfe, 1] + sum(M[ncp, i] * dt * dx[nx, nfe, i] for i in 1:NCP)
        Cont_Eq_First[   nx in 1:Nx],                                     x[nx, 1, 1]       == x0[nx]
        Cont_Eq_rest[    nx in 1:Nx,  nfe in 2:NFE],                      x[nx, nfe, 1]     == x[nx, nfe-1, end]
        Coll_Eq_Quad0[                                 ncp in 1:NCP],     q[1, 1, ncp]      == q0  + sum(M[ncp, i] * dq[1, 1, i] for i in 1:NCP)
        Coll_Eq_Quad[                 nfe in 2:NFE,    ncp in 1:NCP],     q[1, nfe, ncp]    == q[1, nfe-1, NCP]  + sum(M[ncp, i] * dt * dq[1, nfe, i] for i in 1:NCP)
end)
##Objective Function===========================================================
@NLobjective(m1, Min,  q[1, end, end])
##===========================================================
optimize!(m1)
JuMP.termination_status(m1)
JuMP.solve_time(m1)
##===========================================================
Solution = JuMP.value.(u)[:,:]
x       = JuMP.value.(x)[:,:,:]
u_plot  = cat(x0[1], Solution[1:NFE], dims = 1)
x1      = cat(x0[1], x[1,:,:], dims = 1)
x2      = cat(x0[2], x[2,:,:], dims = 1)
t_plot = collect(0:NFE)    
##=========================================================== 
#choose backend for plots
plotlyjs()
p11 = plot(t_plot, x1,         label = "X")
p11 = plot!(t_plot, x2,         label = "X2")

p12 = plot(t_plot, u_plot, xaxis = "Time", label = "U", linetype = :steppost)

fig1 = plot(p11, p12, layout = (2, 1))
