##Initialization
include("Battery_Integrator.jl")
ModelTime = 1;
ISO_sec = 2
delta_sei0 = 1e-10;
P_FR_segment            = NaN*ones(1801)
P_FR_segment[1:600]     .=  20
P_FR_segment[601:1200]  .= -20
P_FR_segment[1201:end]  .= 0
# ab  = convert(Int64, 3600/ISO_sec);
# P_FR_segment            = ones(ab+1)




##Test the Integrator
results = Battery_Model(ModelTime, ISO_sec, delta_sei0, P_FR_segment)


##
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