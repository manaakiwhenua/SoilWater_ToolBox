# =============================================================
#		module: plotSmap
# =============================================================
module plotSmap

   import ..cst, ..hydroStruct, ..kunsat, ..option, ..param, ..path, ..reading, ..wrc
   using Plots.PlotMeasures, LaTeXStrings
   using Suppressor
   using Plots; pgfplotsx

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : PLOT_KUNSAT
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function PLOT_KUNSAT(hydroParam, N_SoilSelect, smap; N_Se= 1000)

         hydroParam2 = hydroStruct.HYDROSTRUCT(N_SoilSelect)

         hydroParam2 = deepcopy(hydroParam)

         Ψ_θΨ_Max = 150000.0 + 10000.0
         Ψ_θΨ_Min = 0.0

         RockFragment = collect(eps(100.0):0.5 / N_Se:0.99)

         N = length(RockFragment)

         Kunsat_Sim  = Array{Float64}(undef, (N))
         Kunsat_PeckWatson = Array{Float64}(undef, (N))

         Plot1 = Plots.plot()
         default(titlefont = (20, "times"), legendfontsize = 18, guidefont = (18, :darkgreen), tickfont = (12, :orange), guide = "x", framestyle = :zerolines, yminorgrid = true)

         for iZ in 1:N_SoilSelect
            for iRF = 1:N
               hydroParam2.θs[iZ] = hydroParam.θs[iZ] * (1.0 - RockFragment[iRF])
               hydroParam2.θsMacMat[iZ] = hydroParam.θsMacMat[iZ] * (1.0 - RockFragment[iRF])
               hydroParam2.θr[iZ] = hydroParam.θr[iZ] * (1.0 - RockFragment[iRF])

               Kunsat_Sim[iRF] = kunsat.θΨ_2_KUNSAT(1.0, iZ, hydroParam2, RockFragment[iRF]; TopsoilSubsoil="Topsoil") * cst.MmS_2_CmH 

               Kunsat_PeckWatson[iRF] =kunsat.θΨ_2_KUNSAT(1.0, iZ, hydroParam, 0.0; TopsoilSubsoil="Topsoil") * (2.0 *  (1.0 - RockFragment[iRF]) /  (2.0 + RockFragment[iRF])) * cst.MmS_2_CmH
            end

            Label= smap.Soilname[iZ] * "_" * string(smap.Depth[iZ])
            Plots.plot!(Plot1, RockFragment, Kunsat_Sim, label=Label * "_Pollacco", palette = :darkrainbow)

            Plots.plot!(Plot1, RockFragment,  Kunsat_PeckWatson, line = (:dot, 4), xlabel=L"RockFragment \ [-]", ylabel=L"K_{s} \ [cm \ h ^{-1}]", label=Label * "PeckWatson", palette = :darkrainbow, legend = :outertopright, yaxis=:log, legendtitle = "LEGEND")
         end

         Plots.savefig(Plot1,"D:\\Temp\\Kunsat.svg")

      end  # function: PLOT_KUNSAT

      
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : HYDROPARAM
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function HYDROPARAM0(Ψ_θΨ, Ψ_KΨ, θ_θΨ, N_θΨ, N_SoilSelect, N_KΨ, K_KΨ, Id_Select, hydro, KunsatModel_Lab; N_Se=1000, smap=[])
         println("  ==  START: Plotting HydroParam  ==")
   
         Flag_OtherData = false
   
         θ_Sim       = Array{Float64}(undef, (N_Se))
         θ_OtherData = Array{Float64}(undef, (N_Se))
         Kunsat_Sim  = Array{Float64}(undef, (N_Se))

         Ψ_θΨ_Min = 0.0

         if Flag_OtherData
            Path = "D:\\Main\\MODELS\\SoilWaterToolBox_Main\\Code_SoilWaterToolbox\\INPUT\\DataSoilHydraulic\\Smap20210226\\Smap20210226_ClappHornberger_Constrained_A_Table_ThetaHK.csv"

            # Structure of the hydroparameters
               hydroData = hydroStruct.HYDROSTRUCT(N_SoilSelect)

            # Populate the values of the parameters
               hydroData, ~ = reading.READ_STRUCT(hydroData, Path)
         end # if Flag_OtherData


         for iZ = param.globalparam.N_iZ_Plot_Start: min(param.globalparam.N_iZ_Plot_End, N_SoilSelect)	
            Ψ_θΨ_Max = maximum(Ψ_θΨ[iZ,N_θΨ[iZ]]) + 100000.0

            Ψ_Sim = expm1.(range(log1p(Ψ_θΨ_Min), stop=log1p(Ψ_θΨ_Max), length=N_Se)) 

            θ_θΨ_Max = hydro.Φ[iZ]

            # Simulated 
               for iΨ = 1:N_Se
                  θ_Sim[iΨ] = wrc.Ψ_2_θDual(Ψ_Sim[iΨ], iZ, hydro)

                  if Flag_OtherData
                     θ_OtherData[iΨ] = wrc.Ψ_2_θDual(Ψ_Sim[iΨ], iZ, hydroData)
                  end

                  Kunsat_Sim[iΨ] = kunsat.Ψ_2_KUNSAT(Ψ_Sim[iΨ], iZ, hydro)
               end # iΨ

            # == Title == 
               Title = smap.Soilname[iZ]  * "_" * string(Int64(floor(smap.Depth[iZ]))) * "_" * string(option.hydro.HydroModel)
  

               Title = Title  * "_" * string(option.hydro.σ_2_Ψm)

            # == Ticks ==
               Ticks = Int64.(Ψ_θΨ[iZ,1:N_θΨ[iZ]])
               
            Plot1 = Plots.plot(layout=(2), bottom_margin = 100px, top_margin = 10px, right_margin = 20px, left_margin = 80px, size=(4000,1300))
            default(titlefont = (20, "times"), legendfontsize = 18, guidefont = (18, :darkgreen), tickfont = (12, :blue), guide = "x", framestyle = :zerolines, yminorgrid = true)

            
            #  == Plot_θ_Ψ  ==
               # Plot_θ_Ψ: Observed
                  X = Ψ_θΨ[iZ,1:N_θΨ[iZ]]
                  Y = θ_θΨ[iZ,1:N_θΨ[iZ]]
                  Label = "Obs"
                  Plot_Θψ = Plots.plot!(Plot1, subplot=1, log1p.(cst.Mm_2_kPa * X) ,Y, seriestype=:scatter, label=Label, color=:red, shape=:square, markersize=4)

               # Plot_θ_Ψ: Simulated
                  X = Ψ_Sim[1:N_Se]
                  Y = θ_Sim[1:N_Se]
                  Label = "Sim"
                  Plot_Θψ = Plots.plot!(Plot1, subplot=1, log1p.(cst.Mm_2_kPa*X), Y, seriestype=:line, label=Label, color=:blue, lw=2)

               # Plotting: Other data
               if Flag_OtherData
                  Y = θ_OtherData[1:N_Se]
                  Label = "Smap_Python"
                  Plot_Θψ= Plots.plot!(Plot1, subplot=1,  log1p.(cst.Mm_2_kPa*X), Y, seriestype=:line, label=Label, color=:green, lw=2)
               end

               # Plot_θ_Ψ: Total porosity point
                  X = zeros(Float64,1)
                  X[1] = 0.0
                  Y = zeros(Float64,1)
                  Y[1] = hydro.Φ[iZ]
                  Label = "\$ \\phi \$"
                  Plot_Θψ = Plots.plot!(Plot1, subplot=1, log1p.(cst.Mm_2_kPa * X) , Y, seriestype=:scatter, label= Label, color=:green, shape=:square, markersize=4) 

               # Plot_θ_Ψ: General attributes
                  xlabel!(L"ln(1 \ + \ \psi) \ [kPa]")
                  ylabel!(L"\theta \ [cm^3 cm^{-3}]")
                  Plot_Θψ= Plots.plot!(Plot1, subplot=1, xlims=( log1p.(cst.Mm_2_kPa * Ψ_θΨ_Min), log1p.( cst.Mm_2_kPa * Ψ_θΨ_Max * 1.1)), ylims =(0.0, θ_θΨ_Max), xticks=(log1p.(cst.Mm_2_kPa * Ψ_θΨ[iZ,1:N_θΨ[iZ]]), Int64.(cst.Mm_2_kPa * Ψ_θΨ[iZ,1:N_θΨ[iZ]])), title=Title)

            #  == Plot_K_Ψ  ==
               if option.hydro.KunsatΨ
               # Plot_K_Ψ: Ks
                  K_Ψ_Max = maximum(K_KΨ[iZ,1:N_KΨ[iZ]] )
                  X = zeros(Float64,1)
                  X[1] = 0.0 
                  Y = zeros(Float64,1)
                  Y[1] = K_Ψ_Max
                  Label = "Ks_Max"
                  Plot_kΘ =Plots.plot!(Plot1, subplot=2, log1p.(X.*cst.Mm_2_kPa), log1p.(Y.*cst.MmS_2_CmH), seriestype=:scatter, label=Label , color= :green, shape= :square, markersize=4)


               if option.hydro.KunsatΨ
                  X = Ψ_KΨ[iZ,1:N_KΨ[iZ]]
                  Y = K_KΨ[iZ,1:N_KΨ[iZ]]
                  Label = "Obs"
                  Plot_kΘ = Plots.plot!(Plot1, subplot=2, log1p.(X.* cst.Mm_2_kPa) , log1p.(Y.*cst.MmS_2_CmH), seriestype=:scatter, label=Label, color= :red, shape= :square, markersize=4)
               end # option.hydro.KunsatΨ

               # Plot_K_Ψ: Sim K_Ψ
                  X = Ψ_Sim
                  Y = Kunsat_Sim 
                  Label = "Sim"
                  Plot_kΘ =Plots.plot!(Plot1, subplot=2, log1p.(X.*cst.Mm_2_kPa) , log1p.(Y.*cst.MmS_2_CmH), seriestype=:line, label=Label, color= :blue, lw=2)

               # General attributes
                  Plots.xlabel!(L"ln( 1 \ + \ \psi ) \ [cm]")
                  Plots.ylabel!(L"ln ( 1 \ + \ K (\psi) )  \ [cm \ h^{-1}]")
                  Plot_kΘ =Plots.plot!(Plot1, subplot=2, xlims = (log1p.(cst.Mm_2_kPa*Ψ_θΨ_Min), log1p.(Ψ_θΨ_Max*cst.Mm_2_kPa)), ylims = (log1p(0.0), log1p(K_Ψ_Max* cst.MmS_2_CmH * 1.1)),  xticks=(log1p.(cst.Mm_2_kPa * Ψ_θΨ[iZ,1:N_θΨ[iZ]]), Int64.(cst.Mm_2_kPa * Ψ_θΨ[iZ,1:N_θΨ[iZ]])))
               
               Plot1 = Plots.plot(Plot1,Plot_Θψ, Plot_kΘ)
            else
               Plot1 = Plots.plot(Plot1, Plot_Θψ)
            end # option.hydro.KunsatΨ 

            Path = path.Plots_θΨK * "Lab_ThetaH_" * Title * ".svg" 
            
            # Path = path.Plots_θΨK * "Lab_ThetaH_" * string(path.Model_Name) * "_" * string(Id_Select[iZ]) * ".svg" 
            # Plots.GRBackend()

            
            @suppress begin
               Plots.savefig(Plot1, Path)
            end
            println("    ~  $(Path) ~")
         end # for iZ

         println("  ==  END: Plotting HydroParam  == \n")		
      return nothing
      end  # function: HYDROPARAM

      # =============================================================
      #		module: makie
      # =============================================================
      module makie
         import ..cst, ..hydroStruct, ..kunsat, ..option, ..param, ..path, ..reading, ..wrc, ...readSmap
         using Makie
         using CairoMakie

         CairoMakie.activate!()
      
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         #		FUNCTION : HYDROPARAM
         # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         function HYDROPARAM(Ψ_θΨ, Ψ_KΨ, θ_θΨ, N_θΨ, N_SoilSelect, N_KΨ, K_KΨ, Id_Select, hydro, KunsatModel_Lab; N_Se=1000, smap=[])
            println("  ==  START: Plotting HydroParam  ==")
      
            Flag_OtherData1 = true
            Flag_OtherData2 = false
      
            # ===================== DATA =====================
            θ_Sim             = fill(0.0,N_Se)
            θ_OtherData       = fill(0.0,N_Se)
            θ_OtherData2      = fill(0.0,N_Se)
            θ_OtherData3      = fill(0.0,N_Se)
            Kunsat_Sim        = fill(0.0,N_Se)
            Kunsat_OtherData  = fill(0.0,N_Se)
            Kunsat_OtherData2 = fill(0.0,N_Se)
            Kunsat_OtherData3 = fill(0.0,N_Se)
            θobs =[]

            Ψ_θΨ_Min = 0.0

            if Flag_OtherData1
               Path = "D:\\Main\\MODELS\\SoilWater-ToolBox2\\src\\INPUT\\DataSoilHydraulic\\Smap20210226\\Smap20210226_ClappHornberger_Constrained_A_Table_ThetaHK.csv"

               option.hydro.HydroModel = :ClappHornberger
                  # Structure of the hydroparameters
                     hydroData = hydroStruct.HYDROSTRUCT(N_SoilSelect)  
                  # Populate the values of the parameters
                  option.hydro.HydroModel = :ClappHornberger
                     hydroData, ~ = reading.READ_STRUCT(hydroData, Path)

               Path = "D:\\Main\\MODELS\\SoilWater-ToolBox2\\src\\INPUT\\DataSoilHydraulic\\Smap20210226\\Smap20210226_Loam.csv"
                  option.hydro.HydroModel = :ClappHornberger
                  # Structure of the hydroparameters
                     hydroData2 = hydroStruct.HYDROSTRUCT(N_SoilSelect)
                  # Populate the values of the parameters
                     hydroData2, ~ = reading.READ_STRUCT(hydroData2, Path) 

                 Path =  "D:\\Main\\MODELS\\SoilWater-ToolBox2\\src\\INPUT\\DataSoilHydraulic\\Smap20210226\\Smap20210226_VangenuchtenJules_Constrained_A_Table_ThetaHK.csv"
                  option.hydro.HydroModel = :VangenuchtenJules
                  # Structure of the hydroparameters
                     hydroData3 = hydroStruct.HYDROSTRUCT(N_SoilSelect)
                  # Populate the values of the parameters
                     hydroData3, ~ = reading.READ_STRUCT(hydroData3, Path) 

            end # if Flag_OtherData

            for iZ = param.globalparam.N_iZ_Plot_Start: min(param.globalparam.N_iZ_Plot_End, N_SoilSelect)	
               Ψ_θΨ_Max = maximum(Ψ_θΨ[iZ,N_θΨ[iZ]]) + 100000.0

               Ψ_Sim = expm1.(range(log1p(Ψ_θΨ_Min), stop=log1p(Ψ_θΨ_Max), length=N_Se)) 

               θ_θΨ_Max = hydro.Φ[iZ]

               # Simulated 
                  for iΨ = 1:N_Se
                     option.hydro.HydroModel = :Vangenuchten
                      θ_Sim[iΨ] = wrc.Ψ_2_θDual(Ψ_Sim[iΨ], iZ, hydro)
                      Kunsat_Sim[iΨ] = kunsat.Ψ_2_KUNSAT(Ψ_Sim[iΨ], iZ, hydro)

                     if Flag_OtherData1
                        # ClappHornberger model Smap_Hydro
                        option.hydro.HydroModel = :ClappHornberger
                           θ_OtherData[iΨ] = wrc.Ψ_2_θDual(Ψ_Sim[iΨ], iZ, hydroData)
                           option.hydro.HydroModel = :ClappHornberger
                           Kunsat_OtherData[iΨ] = kunsat.Ψ_2_KUNSAT(Ψ_Sim[iΨ], iZ, hydroData)

                        # ClappHornberger Loam
                        option.hydro.HydroModel = :ClappHornberger
                           θ_OtherData2[iΨ] = wrc.Ψ_2_θDual(Ψ_Sim[iΨ], iZ, hydroData2)
                            option.hydro.HydroModel = :ClappHornberger
                           Kunsat_OtherData2[iΨ] = kunsat.Ψ_2_KUNSAT(Ψ_Sim[iΨ], iZ, hydroData2)
          
                        # VanGenuchten_Jules
                        option.hydro.HydroModel = :VangenuchtenJules
                         θ_OtherData3[iΨ] = wrc.Ψ_2_θDual(Ψ_Sim[iΨ], iZ, hydroData3)
                         option.hydro.HydroModel = :VangenuchtenJules
                        Kunsat_OtherData3[iΨ] =  kunsat.Ψ_2_KUNSAT(Ψ_Sim[iΨ], iZ, hydroData3)

                         θobs₀ =[ [ 0.456,	0.35,	0.28,	0.16],
                        [0.4465,	0.32,	0.25,	0.15],
                        [0.4465,	0.3,	0.17,	0.1],
                        [0.646,	0.46,	0.36,	0.25],
                        [0.6745,	0.48,	0.36,	0.27],
                        [0.6935,	0.52,	0.43,	0.29]]

                         θobs = θobs₀[iZ,:][1] 
                     else
                         option.hydro.HydroModel = :VangenuchtenJules
                         θ_Sim[iΨ] = wrc.Ψ_2_θDual(Ψ_Sim[iΨ], iZ, hydro)

                         option.hydro.HydroModel = :VangenuchtenJules
                        Kunsat_Sim[iΨ] = kunsat.Ψ_2_KUNSAT(Ψ_Sim[iΨ], iZ, hydro)
                     end
                  end # iΨ

               # == Title == 
                  Title = smap.Soilname[iZ]  * "_" * string(Int64(floor(smap.Depth[iZ]))) * "_" * string(option.hydro.HydroModel)

                  Title = smap.Soilname[iZ]  * "  " * string(Int64(floor(smap.Depth[iZ]))) * " mm"
                  # Title = Title  * "_" * string(option.hydro.σ_2_Ψm)

      
               #  ===================== PLOTTING =====================
                  Fig = Makie.Figure()
                              
               #  == Plot_θ_Ψ  ==
                  # Plot_θ_Ψ: General attributes

                  Fig = Figure(backgroundcolor = RGBf0(0.98, 0.98, 0.98),
                     resolution = (1000, 700))

                  Axis1 = Makie.Axis(Fig[1,1], resolution = (1000, 700))

                  Makie.xlims!(Axis1, log1p.(cst.Mm_2_kPa * Ψ_θΨ_Min), log1p.(cst.Mm_2_kPa * Ψ_θΨ_Max * 1.1))
                  Makie.ylims!(Axis1, 0.0, 0.75)
   
                  Axis1.xticks = (log1p.(cst.Mm_2_kPa * Ψ_θΨ[iZ,1:N_θΨ[iZ]]), string.(Int64.(cst.Mm_2_kPa * Ψ_θΨ[iZ,1:N_θΨ[iZ]])))
                  Axis1.xlabel = "ln(1 + Ψ) [kPa]"
                  Axis1.ylabel =  "θ [mm³ mm⁻³]"
                  Axis1.title = Title
                  Axis1.titlesize= 24

                  P_Smap= Makie.scatter!(Fig[1,1], log1p.(cst.Mm_2_kPa .* Ψ_θΨ[iZ,1:N_θΨ[iZ]]), Float64.(θ_θΨ[iZ,1:N_θΨ[iZ]]), color=:blue, markersize=15, marker = '■', label="Smap")

                  # Plot_θ_Ψ: Simulated
                     # Makie.lines!(Fig[1,1], log1p.(cst.Mm_2_kPa .* Ψ_Sim[1:N_Se]), θ_Sim[1:N_Se], color=:blue, linewidth=2, label="Sim")


                  if Flag_OtherData1
                     P_ClappHonb = Makie.lines!(Fig[1,1], log1p.(cst.Mm_2_kPa .* Ψ_Sim[1:N_Se]), θ_OtherData[1:N_Se], color=:red, linewidth=3)

                     P_ClappHonb_Loan = Makie.lines!(Fig[1,1], log1p.(cst.Mm_2_kPa .* Ψ_Sim[1:N_Se]), θ_OtherData2[1:N_Se], color=:yellow1, linewidth=3)

                     P_vangJules = Makie.lines!(Fig[1,1], log1p.(cst.Mm_2_kPa .* Ψ_Sim[1:N_Se]), θ_OtherData3[1:N_Se], color=:green, linewidth=3)

                     Ψ = [0.0 ,1000.0,10000.0,150000.0]
                     P_Lab = Makie.scatter!(Fig[1,1], log1p.(cst.Mm_2_kPa .* Ψ[1:4]), θobs[1:4], color=:darkviolet, linewidth=2, markersize=15)

                     P_Vang = Makie.lines!(Fig[1,1], log1p.(cst.Mm_2_kPa .* Ψ_Sim[1:N_Se]), θ_Sim[1:N_Se], color=:pink, linewidth=3)
                  else
                      Makie.lines!(Fig[1,1], log1p.(cst.Mm_2_kPa .* Ψ_Sim[1:N_Se]), θ_Sim[1:N_Se], color=:blue, linewidth=2)

                  end

                  # Plot_θ_Ψ: Total porosity point
                     X = zeros(Float64,1)
                     X[1] = 0.0
                     Y = zeros(Float64,1)
                     Y[1] = hydro.Φ[iZ]
                     Label = 
                     P_PtotalPorosity = Makie.scatter!(Fig[1,1], log1p.(cst.Mm_2_kPa .* X), Y, color=:slateblue3, markersize=20, marker ="●")
  
                  # Makie.axislegend()

               # == Plot_K_Ψ  ==
               option.hydro.KunsatΨ = true
               if option.hydro.KunsatΨ
                     Axis2 = Makie.Axis(Fig[1,2])
                     Axis2.xticks = (log1p.(cst.Mm_2_kPa * Ψ_θΨ[iZ,1:N_θΨ[iZ]]), string.(Int64.(cst.Mm_2_kPa * Ψ_θΨ[iZ,1:N_θΨ[iZ]])))
                     Yticks = 1:1:6
                     Axis2.yticks = (Yticks,string.(Yticks))
                     Axis2.xlabel = "ln(1 + Ψ) [kPa]"
                     Axis2.ylabel =  "ln ( 1 + K (Ψ) ) [mm h⁻¹]"

                   Pvang = Makie.lines!(Fig[1,2], log1p.(Ψ_Sim[1:N_Se].*cst.Mm_2_kPa), log1p.(Kunsat_Sim[1:N_Se].*cst.MmS_2_MmH), color=:pink, linewidth=3)

                  if Flag_OtherData1

                     Pclapp = Makie.lines!(Fig[1,2], log1p.(Ψ_Sim[1:N_Se].*cst.Mm_2_kPa), log1p.(Kunsat_OtherData[1:N_Se].*cst.MmS_2_MmH), color=:red, linewidth=3)

                     Ploan = Makie.lines!(Fig[1,2], log1p.(Ψ_Sim[1:N_Se].*cst.Mm_2_kPa), log1p.(Kunsat_OtherData2[1:N_Se].*cst.MmS_2_MmH), color=:yellow1, linewidth=3)

                     P9 = Makie.lines!(Fig[1,2], log1p.(Ψ_Sim[1:N_Se].*cst.Mm_2_kPa), log1p.(Kunsat_OtherData3[1:N_Se].*cst.MmS_2_MmH), color=:green, linewidth=3)

                  else
                      Axis2 = Makie.Axis(Fig[1,2])

                     K_Ψ_Max = maximum(K_KΨ[iZ,1:N_KΨ[iZ]])
                     Makie.xlims!(Axis2, log1p.(cst.Mm_2_kPa*Ψ_θΨ_Min), log1p.(Ψ_θΨ_Max*cst.Mm_2_kPa))
                     Makie.ylims!(Axis2,  (log1p(0.0), log1p(K_Ψ_Max* cst.MmS_2_CmH * 1.1)))
                     Axis2.xlabel = "ln(1 + Ψ) [kPa]"
                     Axis2.ylabel =  "ln ( 1 + K (Ψ) ) [cm h⁻¹]"
                       X = zeros(Float64,1)
                     X[1] = 0.0 
                     Y = zeros(Float64,1)
                     Y[1] = K_Ψ_Max
                     Label = "Ks_Max"
                     Makie.scatter!(Fig[1,2], log1p.(X.* cst.Mm_2_kPa) , log1p.(Y.*cst.MmS_2_CmH) ,Y, color=:yellow, markersize=15, marker = '■', label=Label )

                  # PlotK_Ψ: K(Ψ) obs
                     X = Ψ_KΨ[iZ,1:N_KΨ[iZ]]
                     Y = K_KΨ[iZ,1:N_KΨ[iZ]]
                     Label = "Obs"
                     Makie.scatter!(Fig[1,2], log1p.(X.* cst.Mm_2_kPa) , log1p.(Y.*cst.MmS_2_CmH) ,Y, color=:red, markersize=10, marker = '■', label=Label )

                  # Plot_K_Ψ: K(Ψ) sim
                     X = Ψ_Sim
                     Y = Kunsat_Sim 
                     Label = "Sim"
                     Makie.lines!(Fig[1,2], log1p.(Ψ_Sim.*cst.Mm_2_kPa), log1p.(Kunsat_Sim.*cst.MmS_2_CmH), color=:blue, label=Label)
                     
                  end
                  
               end # option.hydro.KunsatΨ 

               # Fig[3, 1] = Makie.Legend(Fig, Axis1, "PLOTS", orientation=:horizontal)

               # Path = path.Plots_θΨK * "Lab_ThetaH_" * Title * ".svg" 

               leg = Fig[1, end+1] = Makie.Legend(Fig, [P_Smap, P_ClappHonb, P_ClappHonb_Loan, P_vangJules, P_Vang, P_Lab, P_PtotalPorosity], ["Smap", "ClapHornberger", "ClapHornberger_Loam", "VangJules", "Vang","Lab","TotalPorosity"])

               Fig[2, 1:2] = leg
               Makie.trim!(Fig.layout)
               leg.orientation = :horizontal
               Makie.trim!(Fig.layout)
               leg.tellheight = true
               
               Path = path.Plots_θΨK * "Lab_ThetaH_" * string(path.Model_Name) * "_" * string(Id_Select[iZ]) * ".svg" 
      
               Makie.save(Path, Fig)
     
               # Displaying figure in VScode
               if option.globalopt.Plot_Show
                  display(Fig)
               end

               println("    ~  $(Path) ~")
            end # for iZ

            println("  ==  END: Plotting HydroParam  == \n")		
         return nothing
         end  # function: HYDROPARAM

      
            
      end  # module: makie
      # ............................................................
   
end  # module: plotSmap
# ............................................................