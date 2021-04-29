##Start Simulation-----------------------------
include("Battery_OCP.jl")
include("Battery_Integrator.jl")



Timeₛᵢₘ = 6;     #Simulation time in hour
Time_Pred = 2;   #prediction horizon time in hour




Πᶠ            = NaN*ones(Timeₛᵢₘ)
Πᵉ            = NaN*ones(Timeₛᵢₘ)
Πᶠ[1:convert(Int64, Timeₛᵢₘ/3)]     .=  3
Πᶠ[convert(Int64, Timeₛᵢₘ/3) + 1:2*convert(Int64, Timeₛᵢₘ/3)]  .= 8
Πᶠ[2*convert(Int64, Timeₛᵢₘ/3) + 1:end]  .= 5
Πᵉ[1:convert(Int64, Timeₛᵢₘ/3)]     .=  8
Πᵉ[convert(Int64, Timeₛᵢₘ/3) + 1:2*convert(Int64, Timeₛᵢₘ/3)]  .= 2
Πᵉ[2*convert(Int64, Timeₛᵢₘ/3) + 1:end]  .= 5
Πᵉ[1:convert(Int64, Timeₛᵢₘ)]     .=  5
Πᶠ_atach = Πᶠ[end]*ones(Time_Pred)
Πᵉ_atach = Πᵉ[end]*ones(Time_Pred)
Πᶠ = cat(Πᶠ, Πᶠ_atach, dims = 1)
Πᵉ = cat(Πᵉ, Πᵉ_atach, dims = 1)



E0 = 4;
Eₒₚₜ = E0;
global Eₒₚₜ

ISO_sec = 1;
N_ISO  = convert(Int64, 3600/ISO_sec);
α_1 =  -1*ones(Timeₛᵢₘ , N_ISO);
α_2  = reshape(α_1[end,:], (1,N_ISO))
α_atach  = repeat(α_2, outer = [Time_Pred, 1])
α = cat(α_1, α_atach, dims = 1)






ModelTime = 1;
delta_sei0 = 1e-10;


##Simulate!

for t in 1:Timeₛᵢₘ
    global  P_to_Battery, results
    #Optimizing E
    
    fr, Ok, Ek, Pk, profit, t_plot = Battery_MPC(Eₒₚₜ,   Πᶠ[t:t + Time_Pred - 1],
                                                            Πᵉ[t:t + Time_Pred - 1], 
                                                            α[t:t + Time_Pred - 1, :],
                                                            ISO_sec);

    P_to_Battery = Pk[1:N_ISO+1]

    #Simulating the Battery for next hour
    results = Battery_Model(ModelTime, ISO_sec, delta_sei0, P_to_Battery)
    Eₒₚₜ = results.u[end][13]
    
end






##*Plotting everything
plotlyjs()
N_plot = size(results.t)[1]

#[results.u[i][1] for i in 1:N_plot]
p1 = Plots.plot(results.t,  [results.u[i][1]  for i in 1:N_plot], label = "csp_avg");
p1 = Plots.plot!(results.t, [results.u[i][2]  for i in 1:N_plot], label = "csp_s");
p1 = Plots.plot!(results.t, [results.u[i][3]  for i in 1:N_plot], label = "csn_avg");
p1 = Plots.plot!(results.t, [results.u[i][4]  for i in 1:N_plot], label = "csn_s");

p2 = Plots.plot(results.t,  [results.u[i][5]  for i in 1:N_plot], label = "I_int");
p2 = Plots.plot!(results.t, [results.u[i][6]  for i in 1:N_plot], label = "phi_p");
p2 = Plots.plot!(results.t, [results.u[i][7]  for i in 1:N_plot], label = "phi_n");
p2 = Plots.plot!(results.t, [results.u[i][8]  for i in 1:N_plot], label = "V or pot");
p2 = Plots.plot!(results.t, [results.u[i][9]  for i in 1:N_plot], label = "It");
p2 = Plots.plot!(results.t, [results.u[i][10] for i in 1:N_plot], label = "Isei");

p3 = Plots.plot(results.t,  [results.u[i][11] for i in 1:N_plot], label = "delta_sei");
p3 = Plots.plot!(results.t, [results.u[i][15] for i in 1:N_plot], label = "Cf");


p4 = Plots.plot(results.t,  [results.u[i][12] for i in 1:N_plot], label = "cm");
p4 = Plots.plot!(results.t, [results.u[i][13] for i in 1:N_plot], label = "cp");
p4 = Plots.plot!(results.t, [results.u[i][14] for i in 1:N_plot], label = "Q");




##
p1
display(p1)
##
p2

##
p3

##

p4