using JuMP
using Ipopt
using Plots
using Dates
using Gurobi
function Battery_MPC(E0, Πᶠ, Πᵉ, α, ISO_sec)

      ##Required Parameters
  
      Plotting = true;
      Save_Figs = false;
  
  
  
      # Prices for FR and DAM
      #Πᶠ            = NaN*ones(N_hour)
      #Πᵉ            = NaN*ones(N_hour)
  
      # Πᶠ[1:convert(Int64, N_hour/3)]     .=  3
      # Πᶠ[convert(Int64, N_hour/3) + 1:2*convert(Int64, N_hour/3)]  .= 8
      # Πᶠ[2*convert(Int64, N_hour/3) + 1:end]  .= 5
  
      # Πᵉ[1:convert(Int64, N_hour/3)]     .=  8
      # Πᵉ[convert(Int64, N_hour/3) + 1:2*convert(Int64, N_hour/3)]  .= 2
      # Πᵉ[2*convert(Int64, N_hour/3) + 1:end]  .= 5
      # Πᵉ[1:convert(Int64, N_hour)]     .=  5
  
  
      #Starting Values
      F0 = 0;
      O0 = 0;
      l0 = 0;
      P0 = 0;
      # E0 = 2.0;
  
      #ISO data
      #ISO_sec = 10;
      N_ISO = convert(Int64, 3600/ISO_sec);
      # α = rand(N_hour , N_ISO)
      # α =  -1*ones(N_hour , N_ISO);
      N_hour = size(α, 1);
  
  
      dt = ISO_sec
  
  
  
  
      #E bands and band coeffiecient
      Eₘᵢₙ = 0;
      Eₘₐₓ = 10;
      ΔE = Eₘₐₓ - Eₘᵢₙ;
      τₗ = 0.2;
      τᵤ = 0.2;
      μₗ = 0.3;
      μᵤ = 0.3;
  
  
      ##Model Defining
      #m1 = Model(Ipopt.Optimizer)
      m1 = Model(Gurobi.Optimizer)
  
  
      ##Variable Defining----------------------------------
      @variable(m1,  F[1: N_hour])
      @variable(m1,  O[1: N_hour])
      @variable(m1,  l[1: N_hour ])
  
      @variable(m1,  P[1:N_hour, 1:N_ISO])
      @variable(m1,  E[1:N_hour, 1:(N_ISO+1)])
  
      ##Bounds and Starting Points----------------------------------
      for k in 1:N_hour, s in 1:(N_ISO+1)
          #Bounds
          set_lower_bound(F[k], 0)
          set_upper_bound(F[k], 5)
  
  
          set_lower_bound(O[k], 0)
          set_upper_bound(O[k], 5)
  
  
          set_lower_bound(l[k], 0)
          set_upper_bound(l[k], 5)
  
          set_lower_bound(E[k, s], Eₘᵢₙ)
          set_upper_bound(E[k, s], Eₘₐₓ)
          #Start Values
  
          set_start_value(F[k], F0)
          set_start_value(O[k], O0)
          set_start_value(l[k], l0)
          set_start_value(E[k, s], E0)
      end
  
      for k in 1:N_hour, s in 1:(N_ISO)
          set_lower_bound(P[k, s], -5)
          set_upper_bound(P[k, s], 5)
          set_start_value(P[k, s], P0)
      end
  
      ##Constraints----------------------------------
  
      #Initialization of differential state(s)
      @constraints(m1, begin
          Initial_State, E[1, 1] == E0
      end)
      #Continuity
      @constraints(m1, begin
          Connecting[k in 1:N_hour-1], E[k+1, 1] == E[k, end]
      end)
  
      #C1
      @constraints(m1, begin
          Battery_Power[k in 1:N_hour, s in 1:N_ISO], P[k, s] == α[k, s] * F[k] + O[k] - l[k]
      end)
  
      #C2
      @constraints(m1, begin
          Battery_Energy[k in 1:N_hour, s in 1:N_ISO], E[k, s+1] == E[k, s] + P[k,s] * dt/3600
      end)
  
      #C3
      @constraints(m1, begin
          Battery_Safety[k in 1:N_hour, s in 1:(N_ISO+1)], Eₘᵢₙ + τₗ*ΔE <= E[k, s] <= Eₘₐₓ - τᵤ*ΔE
      end)
  
      #C4
      @constraints(m1, begin
          Battery_Terminate[k in 1:N_hour], Eₘᵢₙ + μₗ*ΔE <= E[end, end] <= Eₘₐₓ - μᵤ*ΔE
      end)
  
  
      ##Objective Function----------------------------------
      @objective(m1, Min,  sum( (Πᵉ[k] * O[k]) - (Πᶠ[k] * F[k]) for k in 1:N_hour))
      ##Optimization!----------------------------------
      optimize!(m1)
      JuMP.termination_status(m1)
      JuMP.solve_time(m1)
  
  
      ##Extracting Data----------------------------------
  
      plotlyjs()
      #Creating time axis
      #t_plot = collect(range(0, stop = N_hour, length = N_ISO*N_hour ) ) #!FIX STEPS CORRECTLY
  
      if Plotting
          global t_plot, fr, Ok, Ek, Pk, profit 
          t_plot = Float64[];
          
          for i in 1:N_hour
            if i != N_hour
                  t_sample = collect(range(i-1, stop = i, length = N_ISO+1 ) );
                  pop!(t_sample);
            else
                  t_sample = collect(range(i-1, stop = i, length = N_ISO ) );
            end
            
            t_plot = cat(t_plot, t_sample, dims = 1)
            end
  
  
          fr = JuMP.value.(F)[:]
          Ok = JuMP.value.(O)[:]
          Ek = JuMP.value.(E)[:,1:end-1]
          Pk = JuMP.value.(P)[:,:]
  
          Ek = reshape(transpose(Ek), (1,:))
          Pk = reshape(transpose(Pk), (1,:))
  
          profit = Πᶠ.*fr .- Πᵉ.*Ok;
  
  
          p11 = plot(t_plot, Ek[1,:], label = "E")
          p11 = plot!(t_plot, Pk[1,:], label = "P")
          p11 = plot!( xlims=(0,N_hour+0.5))
          #t_plotₖ = cat(t_plot[N_ISO+1:N_ISO:end], t_plot[end], dims = 1);
          t_plotₖ = t_plot[1:N_ISO:end]
  
  
          p12 = plot(  t_plotₖ, fr, ls = :dot,  markershape = :circle, mc = :match,label = "F")
          p12 = plot!( t_plotₖ, Ok, ls = :dot, markershape = :circle, mc = :match,label = "O")
          p12 = plot!( xlims=(0,N_hour+0.5), ylims=(-1,6))
  
          p13 = plot(  t_plotₖ, Πᶠ, ls = :dot,  markershape = :hexagon, label = "Πᶠ")
          p13 = plot!( t_plotₖ, Πᵉ, ls = :dot,  markershape = :hexagon,label = "Πᵉ")
          #p13 = plot!( t_plotₖ, profit, ls = :dot,  markershape = :hexagon, mc = :match,label = "Profit")
          p13 = plot!( xlims=(0,N_hour+0.5), ylims=(1,9))
  
  
          p14 = plot( t_plotₖ, profit, ls = :dot,  markershape = :hexagon, mc = :match,label = "Profit")
          p14 = plot!( xlims=(0,N_hour+0.5), ylims=(-18,40))
  
  
          p2 = plot(p12, p13, p11, p14, layout = (4, 1), legend=:topright)
  
          #display(p11)
          #display(p2)
  
      end
  
      ##
  
      if Save_Figs
          
          ins = today();
          savefig(p2,"~/Documents/Julia/1.Thesis/Battery_FP/Figures/OCP_Result_$ins.eps")
      end
  
  
  
  
  
  
  
      return fr, Ok, Ek, Pk, profit, t_plot
end