##Initialization

E0 = 4;
##
N_hour = 6;
Πᶠ            = NaN*ones(N_hour)
Πᵉ            = NaN*ones(N_hour)

Πᶠ[1:convert(Int64, N_hour/3)]     .=  3
Πᶠ[convert(Int64, N_hour/3) + 1:2*convert(Int64, N_hour/3)]  .= 8
Πᶠ[2*convert(Int64, N_hour/3) + 1:end]  .= 5

Πᵉ[1:convert(Int64, N_hour/3)]     .=  8
Πᵉ[convert(Int64, N_hour/3) + 1:2*convert(Int64, N_hour/3)]  .= 2
Πᵉ[2*convert(Int64, N_hour/3) + 1:end]  .= 5
Πᵉ[1:convert(Int64, N_hour)]     .=  5



ISO_sec = 10;
N_ISO  = convert(Int64, 3600/ISO_sec);
α =  -1*ones(6 , N_ISO)
##Adding Optimizer
include("Battery_OCP.jl")

##
fr, Ok, Ek, Pk, profit, t_plot = Battery_MPC(E0, Πᶠ, Πᵉ, α, ISO_sec)

