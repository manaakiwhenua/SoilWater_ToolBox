# =============================================================
#		module: plotOther
# =============================================================
module plotOther

		import ..wrc,  ..ΨminΨmax, ..hydroRelation, ..hydroStruct, ..param, ..tool, ..cst, ..kunsat
		using Plots.PlotMeasures, LaTeXStrings
		using Plots; pgfplotsx()
		export ΨMINΨMAX, WOF_STEPS, SE_Ψ_CONSTRAINED, PLOT_σ_2_θr, PLOT_θΨ_Δθ, σ_ψM_SCEARIO


	# ========================================
	# 		ΨMINΨMAX
	# ======================================== 
		function ΨMINΨMAX(hydro, pathHyPix)
			# PREPARING THE DATA
				N_Se = 1000
				local θplot = fill(0.0::Float64, N_Se)
			N_σ = 5
			σ = range(hydro.σ_Min[1], stop=hydro.σ_Max[1], length=N_σ)
         Ψm_Unique = fill(0.0::Float64, N_σ)

         Ψ_Min_σ   = fill(0.0::Float64, N_σ)
         Ψ_Max_σ   = fill(0.0::Float64, N_σ)

			hydroHorizon₂ = hydroStruct.HYDROSTRUCT(N_σ)

			Plots.PGFPlotsXBackend()
			Plot1=Plots.plot(layout=1)

			for iσ = 1:N_σ
				Ψm_Unique[iσ] = hydroRelation.σ_2_Ψm(σ[iσ], param.hydro.kg.Ψσ, param.hydro.kg.Ψσ_Min, param.hydro.kg.Ψσ_Max)
				Ψ_Max_σ[iσ], Ψ_Min_σ[iσ] = ΨminΨmax.ΨMINΨMAX(1.0, 0.85, σ[iσ], param.hydro.kg.σMac, Ψm_Unique[iσ], param.hydro.kg.ΨmMac; Pσ=3.0)
			
				Ψplot = exp.(range(log(Ψ_Min_σ[iσ]), stop=log(Ψ_Max_σ[iσ]), length=N_Se)) 

				hydroHorizon₂.θs[iσ] = 1.0
				hydroHorizon₂.θsMacMat[iσ] = 0.85
				hydroHorizon₂.θr[iσ] = 0.0
				hydroHorizon₂.σ[iσ] = σ[iσ]
				hydroHorizon₂.Ψm[iσ] = Ψm_Unique[iσ]
				hydroHorizon₂.ΨmMac[iσ] = param.hydro.kg.ΨmMac
				hydroHorizon₂.σMac[iσ] = param.hydro.kg.σMac
				
				for iΨ = 1:N_Se
					θplot[iΨ] = wrc.Ψ_2_θDual(Ψplot[iΨ],iσ, hydroHorizon₂)			
				end # for iΨ

				Ψm_Model = wrc.Se_2_ΨDual(0.5, iσ, hydroHorizon₂)

				Plots.plot!(Plot1, log.(Ψplot), θplot, palette=:darkrainbow, label="\$ \\sigma = $(round(σ[iσ],digits=2)) \\ ; ln \\ \\Psi m = $(ceil(Int,log(Ψm_Unique[iσ]))) \$" )
				if iσ == N_σ
					Plots.scatter!(Plot1, log.(Ψ_Min_σ[iσ:iσ]), [[1.0]], marker=(:dtriangle, 7, 1.0,:red), label="\$ \\Psi min \$")
					Plots.scatter!(Plot1, log.([Ψm_Model]), [[0.5]], marker=(:diamond, 7, 1.0,:green), label="\$ \\theta (\\Psi = \\Psi m) \$")

					Plots.scatter!(Plot1, log.(Ψ_Max_σ[iσ:iσ]), [[0.0]],  marker=(:utriangle, 7, 1.0, :blue), label="\$ \\Psi max \$")
				else
					Plots.scatter!(Plot1, log.(Ψ_Min_σ[iσ:iσ]), [[1.0]], label=false, marker=(:dtriangle, 7, 1.0,:red))
					Plots.scatter!(Plot1, log.([Ψm_Model]), [[0.5]], label=false, marker=(:diamond, 7, 1.0,:green))
					Plots.scatter!(Plot1, log.(Ψ_Max_σ[iσ:iσ]), [[0.0]],  label=false, marker=(:utriangle, 7, 1.0, :blue))
				end
			end # iσ ............................................
			# , xlims =(log(minimum(Ψ_Min_σ)), log(maximum(Ψ_Max_σ)))
			Plots.plot!(Plot1,  xlabel=L"ln \ \psi \ [mm]", ylabel=L"Normalised \ \theta \ [mm^3 mm^{-3}]", ylims =(0.0, 1.0), size=(1000,400), legend=:topright, framestyle = [:box :semi :origin :zerolines :grid :true])

			Plots.savefig(Plot1, pathHyPix.Plot_Ψmin_Ψmax)
			println("			 ~ ", pathHyPix.Plot_Ψmin_Ψmax, "~")

		end # function: ΨMINΨMAX



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : WOF_STEPS
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function WOF_STEPS(pathHyPix)

			Path_Output = pathHyPix.Plot_OfStep * "Multiplestep.svg"

			rm(Path_Output, force=true, recursive=true)	

			Label = ["Waitoa";"Otorohanga";"Taupo";"Waihou";"Hamilton"]


			Plot1=Plots.plot(layout=(2,2), size=(1000,600), bottom_margin=20px, grid=:x)

			# TICKS 
				Ticks =[ "1","2a","2b","3a","3b","4a","4b","5a","5b"]

			# OF STEP
				Path = pathHyPix.Input_OfStep * "Of_Step.CSV"

				Of_Waitoa, N_Waitoa = tool.readWrite.READ_HEADER(Path, "Waitoa")
				Of_Waihou, N_Waihou = tool.readWrite.READ_HEADER(Path, "Waihou")
				Of_Taupo, N_Otorohanga = tool.readWrite.READ_HEADER(Path, "Taupo")
				Of_Otorohanga, N_Otorohanga = tool.readWrite.READ_HEADER(Path, "Otorohanga")
				Of_Hamilton, N_Hamilton = tool.readWrite.READ_HEADER(Path, "Hamilton")
				
				Id=1:N_Waitoa

				Of = [Of_Waitoa Of_Waihou Of_Taupo Of_Otorohanga Of_Hamilton]

				for i=1:5
					Plots.plot!(Plot1, subplot=1, Id, Of[1:N_Waitoa, i], palette=:darkrainbow,  marker=(:circle, 4, 1.0), line=(2.0,:solid))
				end
				Plots.plot!(Plot1, subplot=1,  xlabel="", ylabel=L"WOF _{\theta} \ [mm^{3} \ mm^{-3}]", xticks=(1:1:9, Ticks), xtickfont=(12, :white), legend=false, title="(a) Weighted Objective Function", titlelocation = :left)

			# GROUNDWATER STEP
				Path = pathHyPix.Input_OfStep * "Groundwater_Step.csv"

				Groundwater_Waitoa, N_Waitoa = tool.readWrite.READ_HEADER(Path, "Waitoa")
				Groundwater_Waihou, N_Waihou = tool.readWrite.READ_HEADER(Path, "Waihou")
				Groundwater_Taupo, N_Otorohanga = tool.readWrite.READ_HEADER(Path, "Taupo")
				Groundwater_Otorohanga, N_Otorohanga = tool.readWrite.READ_HEADER(Path, "Otorohanga")
				Groundwater_Hamilton, N_Hamilton = tool.readWrite.READ_HEADER(Path, "Hamilton")

				Groundwater = [Groundwater_Waitoa Groundwater_Waihou Groundwater_Taupo Groundwater_Otorohanga Groundwater_Hamilton]

				for i=1:5
					Plots.plot!(Plot1, subplot=2, Id, Groundwater[1:N_Waitoa, i], palette=:darkrainbow,  marker=(:circle, 4, 1.0), line=(2.0,:solid))
				end
				# Plots.plot!(Plot1, subplot=2, xlabel=L"Multistep \ Optimisation \ Steps", ylabel=L"\zeta _{Q} ", xticks=(1:1:8, Ticks), legend=false, title="Groundwater", titlelocation = :left)

				Plots.plot!(Plot1, subplot=2,  xlabel="", ylabel=L"\zeta _{Q} \ [\%]", xticks=(1:1:9, Ticks), xtickfont=(12, :white), legend=false, title="(b) Drainage", titlelocation = :left)


			# EvapoTranspiration STEP
				Path = pathHyPix.Input_OfStep * "Sink_Step.csv"

				Sink_Waitoa, N_Waitoa = tool.readWrite.READ_HEADER(Path, "Waitoa")
				Sink_Waihou, N_Waihou = tool.readWrite.READ_HEADER(Path, "Waihou")
				Sink_Taupo, N_Otorohanga = tool.readWrite.READ_HEADER(Path, "Taupo")
				Sink_Otorohanga, N_Otorohanga = tool.readWrite.READ_HEADER(Path, "Otorohanga")
				Sink_Hamilton, N_Hamilton = tool.readWrite.READ_HEADER(Path, "Hamilton")

				Sink = [Sink_Waitoa Sink_Waihou Sink_Taupo Sink_Otorohanga Sink_Hamilton]

				for i=1:5
					Plots.plot!(Plot1, subplot=3, Id, Sink[1:N_Waitoa, i],  palette=:darkrainbow,  marker=(:circle, 4, 1.0), line=(2.0,:solid))
				end
				# Plots.plot!(Plot1, subplot=3, xlabel=L"Multistep \ Optimisation \ Steps?", ylabel=L"\zeta _{et} ", xticks=(1:1:8, Ticks), legend=false, title="EvaopoTranspiration", titlelocation = :left)

				Plots.plot!(Plot1, subplot=3, xlabel=L"Multistep \ optimisation \ [Layers]", ylabel=L"\zeta _{et}  \ [\%]", xticks=(1:1:9, Ticks), tickfont=(12, :black), legend=false, title="(c) Evapotranspiration", titlelocation = :left)


			# Soil Water Content STEP
				Path = pathHyPix.Input_OfStep * "Swc_Step.csv"

				Swc_Waitoa, N_Waitoa = tool.readWrite.READ_HEADER(Path, "Waitoa")
				Swc_Waihou, N_Waihou = tool.readWrite.READ_HEADER(Path, "Waihou")
				Swc_Taupo, N_Otorohanga = tool.readWrite.READ_HEADER(Path, "Taupo")
				Swc_Otorohanga, N_Otorohanga = tool.readWrite.READ_HEADER(Path, "Otorohanga")
				Swc_Hamilton, N_Hamilton = tool.readWrite.READ_HEADER(Path, "Hamilton")

				Swc = [Swc_Waitoa Swc_Waihou Swc_Taupo Swc_Otorohanga Swc_Hamilton]

				for i=1:5
					Plots.plot!(Plot1, subplot=4, Id, Swc[1:N_Waitoa, i], palette=:darkrainbow, marker=(:circle, 4, 1.0), label=Label[i], line=(2.0,:solid))
				end
				# Plots.plot!(Plot1, subplot=4, xlabel=L"Multistep \ Optimisation \ Steps", ylabel=L"\zeta _{\theta} ", xticks=(1:1:8, Ticks), legend=(-0.15,-0.18), title="Soil Water Content", titlelocation = :left)

				Plots.plot!(Plot1, subplot=4, xlabel=L"Multistep \ optimisation \ [Layers]", ylabel=L"\zeta_{swc}  \ [\%]", xticks=(1:1:9, Ticks), tickfont=(12, :black), legend=(0.75,1.0), title="(d) Root zone soil water content", titlelocation = :left)

				Plots.savefig(Plot1, Path_Output)
				println("			 ~ ", Path_Output, "~")
		end  # function: WGroundwater_STEPS



	# ========================================
	# PLOTTING HYDRAULIC RELATIONSHIP FOR EVERY HORIZON
	# ======================================== 
		function SE_Ψ_CONSTRAINED(hydro, pathHyPix)
				
			# PREPARING THE DATA
				N_Se = 1500
				local θplot_Min    = fill(0.0::Float64, N_Se)
				local θplot_Max    = fill(0.0::Float64, N_Se)
				local Kplot_Min    = fill(0.0::Float64, N_Se)
				local Kplot_Max    = fill(0.0::Float64, N_Se)
		
			N_σ = 3
			σ = range(hydro.σ_Min[1], stop=hydro.σ_Max[1], length=N_σ)

			Ψm_Unique_Min = fill(0.0::Float64, N_σ)
			Ψm_Unique_Max = fill(0.0::Float64, N_σ)
			for iσ = 1:N_σ
				Ψm_Unique_Min[iσ] = hydroRelation.σ_2_Ψm(σ[iσ], param.hydro.kg.Ψσ_Min, param.hydro.kg.Ψσ_Min, param.hydro.kg.Ψσ_Max)
				Ψm_Unique_Max[iσ] = hydroRelation.σ_2_Ψm(σ[iσ], param.hydro.kg.Ψσ_Max, param.hydro.kg.Ψσ_Min, param.hydro.kg.Ψσ_Max)
			end

			# Min and Max Ψ
			Ψ_Min_σ_Min = fill(0.0::Float64, N_σ)
			Ψ_Max_σ_Min = fill(0.0::Float64, N_σ)
			Ψ_Min_σ_Max = fill(0.0::Float64, N_σ)
			Ψ_Max_σ_Max = fill(0.0::Float64, N_σ)
			for iσ=1:N_σ
				Ψ_Min_σ_Min[iσ] , Ψ_Max_σ_Min[iσ] = ΨminΨmax.ΨMINΨMAX(1.0, 0.85, σ[iσ], param.hydro.kg.σMac, Ψm_Unique_Min[iσ], param.hydro.kg.ΨmMac)

				Ψ_Min_σ_Max[iσ] , Ψ_Max_σ_Max[iσ] = ΨminΨmax.ΨMINΨMAX(1.0, 0.85, σ[iσ], param.hydro.kg.σMac, Ψm_Unique_Max[iσ], param.hydro.kg.ΨmMac)
			end  # for iZ=1:N_iZ

			Plots.pgfplotsx()
			Plot1=plot(layout=(1,2), size=(1000,400))

			# FOR EVERY HORIZON
			hydroHorizon₂ = hydroStruct.HYDROSTRUCT(N_σ)
			for iZ = 1:N_σ

				Ψplot_Min = exp.(range(log(Ψ_Min_σ_Min[iZ]), stop=log(Ψ_Max_σ_Min[iZ]), length=N_Se))
				Ψplot_Max = exp.(range(log(Ψ_Min_σ_Max[iZ]), stop=log(Ψ_Max_σ_Max[iZ]), length=N_Se)) 

				Ψplot_Min_Ks = exp.(range(log(Ψ_Min_σ_Min[iZ]), stop=log(Ψ_Max_σ_Min[iZ]), length=N_Se)).\10000.0
				Ψplot_Max_Ks = exp.(range(log(Ψ_Min_σ_Max[iZ]), stop=log(Ψ_Max_σ_Max[iZ]), length=N_Se)) .\10000.0

				hydroHorizon₂.θs[iZ] = 1.0
				hydroHorizon₂.θsMacMat[iZ] = 0.85
				hydroHorizon₂.θr[iZ] = 0.0
				hydroHorizon₂.σ[iZ] = σ[iZ]
				
				hydroHorizon₂.ΨmMac[iZ] = param.hydro.kg.ΨmMac
				hydroHorizon₂.σMac[iZ] = param.hydro.kg.σMac
				hydroHorizon₂.Ks[iZ] = 1.0
				
				for iΨ = 1:N_Se
					hydroHorizon₂.Ψm[iZ] = Ψm_Unique_Min[iZ]
					θplot_Min[iΨ]    = wrc.Ψ_2_θDual(Ψplot_Min[iΨ], iZ, hydroHorizon₂)
					Kplot_Min[iΨ]    = kunsat.Ψ_2_KUNSAT(Ψplot_Min_Ks[iΨ], iZ, hydroHorizon₂)		
					
					hydroHorizon₂.Ψm[iZ] = Ψm_Unique_Max[iZ]
					θplot_Max[iΨ]    = wrc.Ψ_2_θDual(Ψplot_Max[iΨ], iZ, hydroHorizon₂)
					Kplot_Max[iΨ]    = kunsat.Ψ_2_KUNSAT(Ψplot_Max_Ks[iΨ], iZ, hydroHorizon₂)		
				end # for iΨ

			# Plot 1: θplot(Ψplot)
			σ₁ = round(σ[iZ],digits=2)
			Ψm_Unique_Ln = Int(ceil(log(Ψm_Unique_Min[iZ])))

			Colour = [:red, :blue, :green]

			Plots.plot!(Plot1, subplot=1, log.(Ψplot_Min), θplot_Min, ylims =(0.0, 1.0), color=Colour[iZ], label="\$ \\sigma = $σ₁ \\ ln \\ \\Psi m = $Ψm_Unique_Ln \$", line=(2), legend=false, xlims=(0, 30 ))


			Plots.plot!(Plot1, subplot=1, log.(Ψplot_Max), θplot_Max, xlabel=L"ln \ \psi \ [mm]", ylabel=L"S_{e} \ [-]", ylims =(0.0, 1.0), color=Colour[iZ], line=(:dot, 2), guidefont = (16, :black), tickfont = (14, :midnightblue), xlims=(0,30), legend=false, grid=false)

			Ψm_Unique_Ln = Int(ceil(log(Ψm_Unique_Min[iZ])))

			Plots.plot!(Plot1, subplot=2, log1p.(Ψplot_Min_Ks), Kplot_Min, color=Colour[iZ], line=(2), label="\$ \\sigma = $σ₁ \\ ; \\ ln \\ \\Psi m = $Ψm_Unique_Ln \$",)

			Ψm_Unique_Ln = Int(ceil(log(Ψm_Unique_Max[iZ])))

			Plots.plot!(Plot1, subplot=2, log1p.(Ψplot_Max_Ks), Kplot_Max, color=Colour[iZ], ylims =(0.0, 1.0), xlims=(0,10), line=(:dot, 2), xlabel=L"ln(1+\psi) \ [mm]", ylabel=L"K_{Se} \ [-]", label="\$ \\sigma = $σ₁ \\ ; \\ ln \\Psi m = $Ψm_Unique_Ln \$", legend=(0.65,1) , legendfontsize = 14, guidefont = (18, :black), tickfont = (14, :midnightblue), grid=false)

	
			end #iZ ............................................
			# ylims=(0.0,1.0),

			Plots.savefig(Plot1, pathHyPix.Plot_Se_Ψ_Constrained)
			println("			 ~ ",  pathHyPix.Plots_θ∂θ∂Ψ, "~")
		end # function: SE_Ψ_CONSTRAINED



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PLOT_θr
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function PLOT_σ_2_θr(hydro, pathHyPix)
			N_σ = 100
			σ = range(hydro.σ_Min[1], stop=hydro.σ_Max[1], length=N_σ)

			hydro₂=deepcopy(hydro)

			θr = fill(0.0::Float64, N_σ)
			for i=1:N_σ
				hydro₂.σ[1] = σ[i]
				θr[i] = hydroRelation.σ_2_θr(hydro₂, 1)
			end

			Plot1=Plots.plot(layout=1, size=(800,400))

			default(titlefont = (20, "times"), legendfontsize = 18, guidefont = (20, "times", :black), tickfont = (16, :midnightblue), guide = "x", framestyle=:origin)

			Plots.plot!(Plot1, σ, θr, label=false, xlabel=L"\sigma \ [-]", ylabel=L"\theta _r \ [m^3 \ m^{-3}]", linewidth = 1.8, linecolour=:mediumblue, grid=false, ticks=true)

			Plots.savefig(Plot1, pathHyPix.Plot_σ2θr)
			println("			 ~ ", pathHyPix.Plot_σ2θr, "~")

		end # function: PLOT_θr



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PLOT_θψ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function PLOT_θΨ_Δθ(hydro, pathHyPix)
		N_Se = 100

		hydroHorizon₂ = hydroStruct.HYDROSTRUCT(1)
		
		hydroHorizon₂.θs[1] = 0.5
		hydroHorizon₂.θsMacMat[1] = hydroHorizon₂.θs[1]
		hydroHorizon₂.θr[1] = 0.001
		hydroHorizon₂.σ[1] = 2.0
		hydroHorizon₂.Ψm[1] = hydroRelation.σ_2_Ψm(hydroHorizon₂.σ[1], param.hydro.kg.Ψσ,  hydro.Ψm_Min[1],  hydro.Ψm_Max[1])
		hydroHorizon₂.ΨmMac[1] = param.hydro.kg.ΨmMac
		hydroHorizon₂.σMac[1] = param.hydro.kg.σMac

		# Feasible range 
		Ψ_Max, Ψ_Min= ΨminΨmax.ΨMINΨMAX(hydroHorizon₂.θs[1], hydroHorizon₂.θsMacMat[1], hydroHorizon₂.σ[1] , param.hydro.kg.σMac, hydroHorizon₂.Ψm[1], param.hydro.kg.ΨmMac)

		# Plotting the curve
			Ψplot = exp.(range(log(Ψ_Min), stop=log(Ψ_Max), length=N_Se)) 
			θplot = fill(0.0::Float64, N_Se)
			for iΨ = 1:N_Se
				θplot[iΨ] = wrc.Ψ_2_θDual(Ψplot[iΨ], 1, hydroHorizon₂)			
			end # for iΨ

		# Plots.PGFPlotsXBackend()
			Plot1=Plots.plot(layout=1, framestyle = :origin, xlims =(log(Ψ_Min), log(Ψ_Max)), grid=false)

			Plots.plot!(Plot1, log.(Ψplot), θplot, label=false, xlabel=L"ln \ \psi \ [mm]", ylabel=L"\theta \ [mm^3 \ mm^{-3}]")
		# Plotting points on the curve
			N_Δθ = 10
			# Ψ_Δθ = fill(0.0::Float64, N_Δθ)
			θ_Min =  wrc.Ψ_2_θDual(Ψ_Min, 1, hydroHorizon₂) +  eps(10.0)
			θ_Max =  wrc.Ψ_2_θDual(Ψ_Max, 1, hydroHorizon₂) -  eps(10.0)
			Δθ = range(θ_Min, stop=θ_Max , length=N_Δθ)
			for iθ = 1:N_Δθ
				Ψ_Δθ = wrc.θ_2_ΨDual(Δθ[iθ], 1, hydroHorizon₂)

				Plots.scatter!(Plot1, [log(Ψ_Δθ),log(Ψ_Δθ)], [Δθ[iθ],Δθ[iθ]], label=false, line=(:dashdot, 1),  color=:green)
				Plots.plot!(Plot1, [log(Ψ_Δθ), log(Ψ_Δθ)], [0.0, Δθ[iθ]], label=false, line=(:dashdot, 1), color=:black)
				Plots.plot!(Plot1, [log(Ψ_Min), log(Ψ_Δθ)], [Δθ[iθ], Δθ[iθ]], label=false, line=(:dashdot, 1), color=:grey)
			end

			Plots.savefig(Plot1, pathHyPix.Plot_θΨ_Δθ)
			println("			 ~ ", pathHyPix.Plot_θΨ_Δθ, "~")
		return
	end  # function: PLOT_θψ


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PLOT_θψ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function σ_ψM_SCEARIO(pathHyPix)
		σ_F, ~   = tool.readWrite.READ_HEADER(pathHyPix.σ_ψM_Scenario, "σ_F")
		σ_G, ~   = tool.readWrite.READ_HEADER(pathHyPix.σ_ψM_Scenario, "σ_G")
		σ_J, ~   = tool.readWrite.READ_HEADER(pathHyPix.σ_ψM_Scenario, "σ_J")
		Ψm_F, ~   = tool.readWrite.READ_HEADER(pathHyPix.σ_ψM_Scenario, "Ψm_F")
		Ψm_G, ~   = tool.readWrite.READ_HEADER(pathHyPix.σ_ψM_Scenario, "Ψm_G")
		Ψm_J, ~   = tool.readWrite.READ_HEADER(pathHyPix.σ_ψM_Scenario, "Ψm_J")

		# Plots.PGFPlotsXBackend()
		Plots.pgfplotsx()
		Plot1=plot(layout=(1,2), size=(1000,400), legend=(0.72,0.15), legendfontsize=14, guidefont=(18, :black), tickfont=(14, :midnightblue), framestyle=:box, grid=false)

		Plots.scatter!(Plot1, subplot=1, σ_F, σ_J, label="Scenario_H", marker=(:circle, 4, 1.0,:blue))
		Plots.scatter!(Plot1, subplot=1, σ_F, σ_G, label="Scenario_K", marker=(:dtriangle, 4, 1.0,:red))
		Plots.plot!(Plot1, subplot=1, [(0,0), (4,4)], label=false, line=(:dashdot, 2), color=:grey)
		Plots.plot!(Plot1, subplot=1, xlabel=L"\sigma \ [-] \ \ Scenario \_ G", ylabel=L"\sigma \ [-]", xlims =(1,4.0), ylims =(1,4.0))

		Plots.scatter!(Plot1, subplot=2, log.(Ψm_F), log.(Ψm_J), label="Scenario_H", marker=(:circle, 4, 1.0,:blue))
		Plots.scatter!(Plot1, subplot=2, log.(Ψm_F), log.(Ψm_G), label="Scenario_K", marker=(:dtriangle, 4, 1.0,:red))
		Plots.plot!(Plot1, subplot=2, [(0,0), (13,13)], label=false, line=(:dashdot, 1), color=:grey)
		Plots.plot!(Plot1, subplot=2 , xlims=(7.0,12.0), ylims =(7.0,12.0))
		Plots.plot!(Plot1, subplot=2, xlims=(7.0,12.0), ylims =(7.0,12.0) , xlabel=L"ln \ \psi_{m} \ [mm] \ Scenario \_ G", ylabel=L"ln \ \psi_{m} \ [mm]")

		Path = pathHyPix.Plots_σΨm * "Plot_σ_ψM_Scenario.svg"
		Plots.savefig(Plot1, Path)
		println("			 ~ ", Path, "~")
	end

end  # module: plotOther
# ............................................................