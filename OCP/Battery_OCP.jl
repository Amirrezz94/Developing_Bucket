using JuMP
using Ipopt
using Plots
##Required Parameters


N_hour = 6;


# Prices for FR and DAM
Πᶠ            = NaN*ones(N_hour)
Πᵉ            = NaN*ones(N_hour)

Πᶠ[1:convert(Int64, N_hour/3)]     .=  0.1
Πᶠ[convert(Int64, N_hour/3) + 1:2*convert(Int64, N_hour/3)]  .= 4
Πᶠ[2*convert(Int64, N_hour/3) + 1:end]  .= 2

Πᵉ[1:convert(Int64, N_hour/3)]     .=  5
Πᵉ[convert(Int64, N_hour/3) + 1:2*convert(Int64, N_hour/3)]  .= 0.1
Πᵉ[2*convert(Int64, N_hour/3) + 1:end]  .= 4


#Starting Values
F0 = 0;
O0 = 0;
l0 = 0;
P0 = 0;
E0 = 0;

#ISO data
ISO_sec = 10;
N_ISO = convert(Int64, 3600/ISO_sec)
α = rand(N_hour , N_ISO)

dt = 1




#E bands and band coeffiecient
Eₘᵢₙ = 0;
Eₘₐₓ = 2;
ΔE = Eₘₐₓ - Eₘᵢₙ;
τₗ = 0.2;
τᵤ = 0.2;
μₗ = 0.3;
μᵤ = 0.3;


##Model Defining
m1 = Model(Ipopt.Optimizer)


##Variable Defining----------------------------------
@variable(m1,  F[1:N_hour])
@variable(m1,  O[1:N_hour])
@variable(m1,  l[1:N_hour])
@variable(m1,  P[1:N_hour, 1:N_ISO])
@variable(m1,  E[1:N_hour, 1:N_ISO])

##Bounds and Starting Points----------------------------------
for k in 1:N_hour, s in 1:N_ISO
    #Bounds
    set_lower_bound(F[k], 0)
    set_upper_bound(F[k], 2)

    set_lower_bound(O[k], 0)
    set_upper_bound(O[k], 2)

    set_lower_bound(P[k, s], 0)
    set_upper_bound(P[k, s], 2)

    set_lower_bound(l[k], 0)
    set_upper_bound(l[k], 2)

    set_lower_bound(E[k, s], Eₘᵢₙ)
    set_upper_bound(E[k, s], Eₘₐₓ)
    #Start Values
    set_start_value(F[k], F0)
    set_start_value(O[k], O0)
    set_start_value(l[k], l0)
    set_start_value(P[k, s], P0)
    set_start_value(E[k, s], E0)
end


##Constraints----------------------------------


#C1
@constraints(m1, begin
      Battery_Power[k in 1:N_hour, s in 1:N_ISO], P[k, s] == α[k, s] * F[k] + O[k] - l[k]
end)


#C2
@constraints(m1, begin
      Battery_Energy[k in 1:N_hour, s in 1:N_ISO-1], E[k, s+1] == E[k, s] + P[k,s] * dt
end)


#C3
@constraints(m1, begin
      Battery_Safety[k in 1:N_hour, s in 1:N_ISO-1], Eₘᵢₙ + τₗ*ΔE <= E[k, s] <= Eₘₐₓ - τᵤ*ΔE
end)
#C4

@constraints(m1, begin
      Battery_Terminate[k in 1:N_hour], Eₘᵢₙ + μₗ*ΔE <= E[k, end] <= Eₘₐₓ - μᵤ*ΔE
end)


##Objective Function----------------------------------
@objective(m1, Min,  sum(- (Πᶠ[k] * F[k]) + (Πᵉ[k] * O[k]) for k in 1:N_hour))
##Optimization!----------------------------------
optimize!(m1)
JuMP.termination_status(m1)
JuMP.solve_time(m1)


##Extracting Data----------------------------------

t_plot = collect(0:dt:dt * N_hour) 
fr = JuMP.value.(F)[:]
Ok = JuMP.value.(O)[:]
Ek = JuMP.value.(E)[:,:]

FR  =   cat(F0, fr[1:N_hour], dims = 1)
DAM_P = cat(O0, Ok[1:N_hour], dims = 1)
E_Plot = cat(E0, Ek[:,end], dims = 1)
Πᶠ_Plot = cat(Πᶠ[1], Πᶠ[1:end], dims = 1)
Πᵉ_Plot = cat(Πᵉ[1], Πᵉ[1:end], dims = 1)
   
##Plotting---------------------------------- 
#choose backend for plots
plotlyjs()

p1 =  plot(t_plot, FR,     label = "Fₖ",    linetype = :steppost);
p1 = plot!(t_plot, DAM_P,  label = "Oₖ", linetype = :steppost);



p2 = plot(t_plot, Πᶠ_Plot, label = "Πᶠ", linetype = :steppost, linestyle = :dash);
p2 = plot!(t_plot, Πᵉ_Plot, label = "Πᵉ", linetype = :steppost, linestyle = :dash);

p3 = plot(t_plot, E_Plot, label = "E", linetype = :steppost, linestyle = :dash);

fig1 = plot(p1, p2, p3, layout = (3, 1), xaxis = "Time")