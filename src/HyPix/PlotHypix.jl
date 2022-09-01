# =============================================================
#		module: plotHypix
# =============================================================
module plotHypix

	import  ..cst, ..kunsat, ..option, ..param, ..rootwateruptake, ..tool, ..wrc, ..ΨminΨmax
	import Dates: value, DateTime
	using PGFPlots

	export θΨK

	

	# ========================================
	# PLOTTING HYDRAULIC RELATIONSHIP FOR EVERY HORIZON
	# ======================================== 
	function θΨK(hydroHorizon, N_iHorizon, iSim, pathHyPix)

		# Deriving the Min and Max Ψ from principals of soil physics
		Ψ_Min_Horizon = fill(0.0::Float64, N_iHorizon)
		Ψ_Max_Horizon = fill(0.0::Float64, N_iHorizon)
		for iZ=1:N_iHorizon
			Ψ_Max_Horizon[iZ], Ψ_Min_Horizon[iZ] = ΨminΨmax.ΨMINΨMAX(hydroHorizon.θs[iZ], hydroHorizon.θsMacMat[iZ], hydroHorizon.σ[iZ], hydroHorizon.σMac[iZ], hydroHorizon.Ψm[iZ], hydroHorizon.ΨmMac[iZ])
		end  # for iZ=1:N_iHorizon
		
		# PREPARING THE DATA
			N_Se = 1000
			local Ψplot = exp.(range(log(minimum(Ψ_Min_Horizon[1:N_iHorizon])), stop = log(maximum(Ψ_Max_Horizon[1:N_iHorizon])), length=N_Se)) 

			local θplot    = fill(0.0::Float64, N_Se)
			local Kplot    = fill(0.0::Float64, N_Se)
			local ∂θ∂Ψplot = fill(0.0::Float64, N_Se)
			local ∂K∂Ψplot = fill(0.0::Float64, N_Se)

			Plot_θΨK = PGFPlots.GroupPlot(4, 100, groupStyle = "horizontal sep = 3.5cm, vertical sep = 3.5cm")

		# FOR EVERY HORIZON
		for iZ = 1:N_iHorizon
			
			for iΨ = 1:N_Se
				if Ψ_Max_Horizon[iZ] ≥ Ψplot[iΨ] ≥ Ψ_Min_Horizon[iZ]
					θplot[iΨ]    = wrc.Ψ_2_θDual(Ψplot[iΨ], iZ, hydroHorizon)
					
					Kplot[iΨ]    = kunsat.Ψ_2_KUNSAT(Ψplot[iΨ], iZ, hydroHorizon)
					
					∂θ∂Ψplot[iΨ] = wrc.∂θ∂Ψ(Ψplot[iΨ], iZ, hydroHorizon)

					∂K∂Ψplot[iΨ] = kunsat.∂K∂Ψ(Ψplot[iΨ], iZ, hydroHorizon)
				else
					θplot[iΨ]    = NaN
					
					Kplot[iΨ]    = NaN
					
					∂θ∂Ψplot[iΨ] = NaN

					∂K∂Ψplot[iΨ] = NaN
				end
			end # for iΨ

			Θs_Max = maximum(hydroHorizon.θs[1:N_iHorizon]) + 0.05
			Ks_Min = 10.0 ^ -7 * cst.MmS_2_CmH
			Ks_Max = maximum(hydroHorizon.Ks[1:N_iHorizon]) * cst.MmS_2_CmH * 1.1

			Title =" $(pathHyPix.SiteName_Hypix)  Layer = $(iZ)"
		
		# Plot 1: θΨ
			Plot_θΨ = PGFPlots.Plots.Linear(log.(Ψplot) , θplot, style=" smooth, blue, very thick", mark="none", legendentry=L"$ \theta ( \Psi ) $")

			Plot_hydro = [Plot_θΨ]

			push!(Plot_θΨK, PGFPlots.Axis(Plot_hydro, style="width=20cm, height=10cm", title=Title, xlabel=L"$ Ln \ \Psi [mm]$", ylabel=L"$ \theta \ [mm{^3} \ mm^{-3}]$", xmin=log(Ψ_Min_Horizon[iZ]), xmax=log(Ψ_Max_Horizon[iZ]), ymin=0.0, ymax=Θs_Max, legendStyle ="{at={(0.0,-0.25)}, anchor=south west, legend columns=1}"))

		# Plot 2: Kplot(Ψplot)
			Plot_Kθ = PGFPlots.Plots.Linear(log.(Ψplot), Kplot .* cst.MmS_2_CmH, style=" smooth, red, very thick", mark="none", legendentry=L"$ K_{unsat} \ ( \Psi ) $")

			Plot_hydro = [Plot_Kθ]
			
			push!(Plot_θΨK, PGFPlots.Axis(Plot_hydro, style="width=20cm, height=10cm", title=Title,  xmin=log(Ψ_Min_Horizon[iZ]), xmax=log(Ψ_Max_Horizon[iZ]), ymax=Ks_Max, ymode="log", xlabel=L"$Ln \  \Psi [mm]$", ylabel=L"$ K_{unsat} \ [cm \ h^{-1}]$", legendStyle ="{at={(0.0,-0.25)}, anchor=south west, legend columns=1}"))

		# Plot 3: ∂θ∂Ψplot
			Plot_∂θ∂Ψ = PGFPlots.Plots.Linear(log.(Ψplot), ∂θ∂Ψplot , style=" smooth, green, very thick", mark="none", legendentry=L"$ \partial \theta \partial \Psi $")

			Plot_hydro = [Plot_∂θ∂Ψ]

			push!(Plot_θΨK, PGFPlots.Axis(Plot_hydro, style="width=20cm, height=10cm", title=Title, xmin=log(Ψ_Min_Horizon[iZ]), xmax=log(Ψ_Max_Horizon[iZ]), xlabel=L"$Ln \  \Psi [mm] $", ylabel=L"$ \partial \theta \partial \Psi $", legendStyle ="{at={(0.0,-0.25)}, anchor=south west, legend columns=1}"))

		# Plot 4: ∂K∂Ψplot
			Plot_∂K∂Ψ = PGFPlots.Plots.Linear(log.(Ψplot), ∂K∂Ψplot, style=" smooth, teal, very thick", mark="none", legendentry=L"$ \partial K \partial \Psi $")

			Plot_hydro = [Plot_∂K∂Ψ]

			push!(Plot_θΨK, PGFPlots.Axis(Plot_hydro, style="width=20cm, height=10cm", title=Title, xmin=log(Ψ_Min_Horizon[iZ]), xmax=log(Ψ_Max_Horizon[iZ]), xlabel=L"$Ln \  \Psi \ [mm]$", ylabel=L"$ \partial K \partial \Psi $", legendStyle ="{at={(0.0,-0.25)}, anchor=south west, legend columns=1}"))

		end #iZ ............................................

		Path = pathHyPix.Plot_Hypix_θΨK * "_" * string(iSim) * ".svg"
		PGFPlots.save(Path, Plot_θΨK) 
	end # function θΨK


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ROOTDENSITY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function VEG_FUNCTIONS(discret, iSim, N_iRoot, veg, Z, ΔRootDensity, pathHyPix)

			Plot_All = PGFPlots.GroupPlot(2, 1, groupStyle = "horizontal sep = 3cm, vertical sep = 3cm")

			# PLOT VEG_FUNCTIONS
				ΔRootDensity_Norm = fill(0.0::Float64, N_iRoot)
				# Taking accoung the tickness of the discretisation
				# for iZ=1:N_iRoot
				# 		ΔRootDensity_Norm[iZ] = Z[N_iRoot] * ΔRootDensity[iZ] / discret.ΔZ[iZ]
				# end

				# Plotting
					Plot_RootDensity = PGFPlots.Plots.Linear(ΔRootDensity[1:N_iRoot], discret.Znode[1:N_iRoot], style=" smooth, cyan, very thick", mark="none")

					Plot = [Plot_RootDensity]

					push!(Plot_All, PGFPlots.Axis(Plot, style="width=12cm, height=8cm", xlabel=L"$ \Delta Rdf \ [\%] $", ylabel=L"$Z \ [mm]$", title="(a)"))
			
			# PLOT StressReduction
				# Data	
				N_Se = 6
				Ψstress = fill(0.0::Float64, 2, N_Se) 
				Ψstress[1,1] = veg.Ψfeddes1 / 10.0
				Ψstress[1,2] = veg.Ψfeddes1
				Ψstress[1,3] = veg.Ψfeddes2
				Ψstress[1,4] = veg.Ψfeddes3
				Ψstress[1,5] = veg.Ψfeddes4
				Ψstress[1,6] = veg.Ψfeddes4 * 2.0

				Wsf = fill(0.0::Float64, N_Se)
				for iΨ ∈ 1:N_Se
					Wsf[iΨ] = rootwateruptake.stressReduction.WATER_STRESS_FUNCTION(2, iΨ, veg, Ψstress)
				end

				Plot_Wsf = PGFPlots.Plots.Linear(Ψstress[1,1:N_Se] .* cst.Mm_2_kPa, Wsf[1:N_Se], style="violet, very thick", mark="none")

				Plot = [Plot_Wsf]

				push!(Plot_All, PGFPlots.Axis(Plot, style="width=12cm, height=8cm", xlabel=L"$ \Psi \ [kPa]$", xmode="log", ylabel=L"$ F_{waterStress} \ [-]$", title="(b)"))

			Path = pathHyPix.Vegetation * "_" * string(iSim) * ".svg"
			PGFPlots.save(Path, Plot_All)	
		end  # function ROOTDENSITY




			# =============================================================
			#		module: plots
			# =============================================================
			module plots
			import ...sorptivity, ..wrc, ..cst, ...option, ..param, ...readHypix
			export PLOT_SORPTIVITY

				using Plots.PlotMeasures, LaTeXStrings
				using Plots;pgfplotsx()

				using Dates
				
				# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				#		FUNCTION : PLOT_SORPTIVITY
				# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				function PLOT_SORPTIVITY(iSim, hydro, pathHyPix)
					println("  ==  START: PLOT_SORPTIVITY_SeIni  ==")

					# Setting the range of values for Se
                  Se_Ini         = collect(0.0:0.001:1.0)
                  N_SeIni        = length(Se_Ini)
                  Sorptivity_Mod = fill(0.0::Float64, (N_SeIni))
                  θ_Ini          = fill(0.0::Float64, (N_SeIni))

					for iSeIni=1:N_SeIni
						θ_Ini[iSeIni] = wrc.Se_2_θ(Se_Ini[iSeIni], 1, hydro)

						Sorptivity_Mod[iSeIni] = sorptivity.SORPTIVITY(θ_Ini[iSeIni], 1, hydro) 
					end
					
					# PLOTTING ====================	
						Plot1=Plots.plot(layout=1)

						Title =" $(pathHyPix.SiteName_Hypix)"

						Plots.plot!(Plot1, Se_Ini[1:N_SeIni] , Sorptivity_Mod[1:N_SeIni], framestyle = [:box :semi :origin :zerolines :grid :true], xlabel=L"Initial \ Se \ [-]", ylabel=L"Sorptivity \  [ \ mm \ \sqrt s \ ]", label="", grid=false) 
					
						Path =pathHyPix.Plot_Sorptivity  * "_" * string(iSim) * ".svg"

						Plots.savefig(Plot1, Path)

						println("			 ~ ", Path, "~")

				end  # function: PLOT_SORPTIVITY


			# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# # 		FUNCTION : INTERCEPTION
			# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# 	function RAINFALL_INTERCEPTION(clim, i∑T_CalibrStart_Day, iSim, pathHyPix)

			# 		# TICKS
			# 			DateTick=range(clim.Date[i∑T_CalibrStart_Day],step=Day(7),clim.Date[clim.N_Climate])
						
			# 			DateTick2= Dates.format.(DateTick, "d u Y")
					
			# 		# PLOT
			# 			Plot1=Plots.plot(layout=1)

			# 			Title =" $(pathHyPix.SiteName_Hypix)" 
						
			# 			Plots.plot!(Plot1, clim.Date[i∑T_CalibrStart_Day:clim.N_Climate], clim.Pr[i∑T_CalibrStart_Day:clim.N_Climate], color=:blue, colorbar=false,  line = :solid, label= L"$\Delta Pr  $")
						
			# 			Plots.plot!(Plot1, clim.Date[i∑T_CalibrStart_Day:clim.N_Climate], clim.Pr_Through[i∑T_CalibrStart_Day:clim.N_Climate], color=:cyan, colorbar=false, label=L"$\Delta Pr_{through}$")

			# 			Plots.plot!(Plot1, clim.Date[i∑T_CalibrStart_Day:clim.N_Climate], 10.0*clim.Pet[i∑T_CalibrStart_Day:clim.N_Climate], color=:green, colorbar=false, label= L"$10x\Delta Pet_{int}$")
						
			# 			Plots.plot!(Plot1, grid=false, framestyle=:origin, size=(1000, 600), legend=:topright, xrotation=rad2deg(pi/3), xticks=(DateTick, DateTick2), title=Title, xlabel=L"$Day$", ylabel=L"$Daily \ \Delta Pr  \ \slash \ \Delta Pr_{through} \ \slash \ \Delta Pet_{int} \ [mm] $")

												
			# 		Path = pathHyPix.Plot_RainfallInterception * "_" * string(iSim) * ".svg"
			# 		Plots.savefig(Plot1, Path)
			# 		println("			 ~ ", Path, "~")
			# 	end  # function: INTERCEPTION

			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#		FUNCTION : TIMESERIES
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				function TIMESERIES(∑T_Date_Plot, ∑T_Plot, obsθ, discret, Flag_Plot_Pond, iSim, N_∑T_Plot, N_iZ, ΔEvaporation_Plot, ΔFlux_Plot, ΔPet_Plot, ΔPond_Plot, ΔPr_Plot, ΔSink_Plot, θ_Plot, θobs_Plot, clim, i∑T_CalibrStart_Day, θsim_Aver, pathHyPix)

				# PATH
					Path = pathHyPix.Plot_HypixTime * "_" * string(iSim) * ".svg"
					rm(Path, force=true, recursive=true)
					
				# READING DATES
					param = readHypix.DATES(pathHyPix)

				# TICKS
					# Date_Start_Calibr = obsθ.Date[1]
					Date_Start_Calibr = DateTime(param.hyPix.obsθ.Year_Start, param.hyPix.obsθ.Month_Start, param.hyPix.obsθ.Day_Start, param.hyPix.obsθ.Hour_Start, param.hyPix.obsθ.Minute_Start, param.hyPix.obsθ.Second_Start) # since we need to compute the culmulativeof the 1rst day
					
					# Date_End_Calibr = obsθ.Date[end]
					Date_End_Calibr = DateTime(param.hyPix.Year_End, param.hyPix.Month_End, param.hyPix.Day_End, param.hyPix.Hour_End, param.hyPix.Minute_End, param.hyPix.Second_End)
					
					DateTick=range(Date_Start_Calibr,step=Day(61),Date_End_Calibr)
					
					DateTick2= Dates.format.(DateTick, "d u Y")
				
				# PLOTTING
					Plot = Plots.plot(layout=(3, 1), size=(2500,2200), bottom_margin=0.01mm)
					
					default(titlefont=(20,"times"), legendfontsize=24, guidefont=18, tickfont=18, grid=true)

				# Plot Climate	
				iSubplot = 0
				if option.hyPix.Plot_Climate
					iSubplot += 1		

					Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, ∑T_Date_Plot[1:N_∑T_Plot], -ΔFlux_Plot[1:N_∑T_Plot, N_iZ+1], label=L"$\Delta Q$", line=(:solid, 1), linecolour=:red, fillcolor=:darkred, fill=(0,:darkred))
					
					Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, ∑T_Date_Plot[1:N_∑T_Plot], ΔPond_Plot[1:N_∑T_Plot], label=L"$\Delta H_{Pond}$", linecolour=:grey, fill = (0, :grey))
					
					Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, clim.Date[i∑T_CalibrStart_Day:clim.N_Climate], clim.Pr[i∑T_CalibrStart_Day:clim.N_Climate], color=:blue, colorbar=false,  line =(:sticks, :solid, 5), label= L"$\Delta Pr  $")

					Plot_Climate = Plots.plot!(Plot, clim.Date[i∑T_CalibrStart_Day:clim.N_Climate], clim.Pr_Through[i∑T_CalibrStart_Day:clim.N_Climate], color=:cyan, line =(:sticks, :solid, 4), colorbar=false, label=L"$\Delta Pr_{through}$")
	
					Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, ylabel=L"$Daily \ Simulation \ [mm]$", title=pathHyPix.SiteName_Hypix, xtickfont = (0.01, :white), xrotation=rad2deg(pi/2))
				end # if: option.hyPix.Plot_Climate

				# PLOT EVAPOYTRANSPIRATION
					iSubplot += 1	

					Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, ∑T_Date_Plot[2:N_∑T_Plot], ΔPet_Plot[2:N_∑T_Plot], linecolour=:darkgreen, label=L"$\Delta Pet$", line=(2.5,:solid))

					Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, ∑T_Date_Plot[2:N_∑T_Plot], ΔSink_Plot[2:N_∑T_Plot], linecolour=:red, line=(2.0,:solid), label=L"$\Delta Sink$")

					Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, ∑T_Date_Plot[2:N_∑T_Plot], (ΔSink_Plot[2:N_∑T_Plot].-ΔEvaporation_Plot[2:N_∑T_Plot]), label=L"$\Delta Rwu$", linecolour=:blue, line=(2.0,:solid))
					
					Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, ∑T_Date_Plot[2:N_∑T_Plot], ΔEvaporation_Plot[2:N_∑T_Plot], label=L"$\Delta Evap$", linecolour=:purple4, line=(2.0,:solid))

					Plot_Climate = Plots.plot!(Plot, subplot=iSubplot, ylabel=L"$Daily \ Simulation \ [mm]$", xtickfont = (0.01, :white), xrotation=rad2deg(pi/2))

				# PLOT Θ
				if option.hyPix.Plot_θ
					iSubplot += 1

					Style_Hypix = [:red, :darkviolet, :orange, :teal, :blue]

					# Observation θplot obs
					for ithetaObs = 1:obsθ.Ndepth
						# lABEL
							Label_Obs = "Obs=" * string(Int(floor(obsθ.Z[ithetaObs]))) * "mm"

							Label_Sim = "Sim=" * string( Int(floor((discret.Znode[obsθ.ithetaObs[ithetaObs]])))) * "mm"

						# Plotting
							# Plot_θ = Plots.plot!(Plot, subplot=iSubplot, ∑T_Date_Plot[1:N_∑T_Plot], θobs_Plot[1:N_∑T_Plot, iZobs].+param.hypix.calibr.θobs_Uncert, line=(0.5,:solid), linecolour=Style_Hypix[iZobs], label=false)
		
							# Plot_θ = Plots.plot!(Plot, subplot=iSubplot, ∑T_Date_Plot[1:N_∑T_Plot], max.(θobs_Plot[1:N_∑T_Plot, iZobs].-param.hypix.calibr.θobs_Uncert, 0.0), line=(0.5,:solid), linecolour=Style_Hypix[iZobs], label=false)
							
							if option.hyPix.θobs_Average
								Plot_θ = Plots.plot!(Plot, subplot=iSubplot, ∑T_Date_Plot[1:N_∑T_Plot], θobs_Plot[1:N_∑T_Plot, ithetaObs], line=(2.5,:solid), linecolour=Style_Hypix[ithetaObs], label="Obs θaver [0-40cm]")

								Plot_θ = Plots.plot!(Plot, subplot=iSubplot, ∑T_Date_Plot[1:N_∑T_Plot], θsim_Aver[1:N_∑T_Plot], label="Sim θaver [0-40cm]", line=(2.5,:solid), linecolour=:blue)

								Plot_θ = Plots.plot!(Plot, subplot=iSubplot, ∑T_Date_Plot[1:N_∑T_Plot], θ_Plot[1:N_∑T_Plot,4], label="Sim θ=10cm", line=(2.5,:dashdot), linecolour=:darkblue)

								Plot_θ = Plots.plot!(Plot, subplot=iSubplot, ∑T_Date_Plot[1:N_∑T_Plot], θ_Plot[1:N_∑T_Plot,14], label="Sim θ=35cm", line=(2.5,:dashdot), linecolour=:darkblue)
						else
							Plot_θ = Plots.plot!(Plot, subplot=iSubplot, ∑T_Date_Plot[1:N_∑T_Plot], θobs_Plot[1:N_∑T_Plot, iZobs], line=(2.5,:solid), linecolour=Style_Hypix[iZobs], label=Label_Obs)

							Plot_θ = Plots.plot!(Plot, subplot=iSubplot, ∑T_Date_Plot[1:N_∑T_Plot], θ_Plot[1:N_∑T_Plot, calibr.iZobs[iZobs]], label=Label_Sim, line=(2.5,:dashdot), linecolour=Style_Hypix[iZobs])
						end  # if: option.hyPix.

					end # loop

					Plot_θ = Plots.plot!(subplot=iSubplot, ylabel=L"$\theta \ [mm^3 \ mm^{-3}]$")

					Plot = Plots.plot(Plot, Plot_θ, Plot_Climate, xmin=∑T_Date_Plot[1], xmax=∑T_Date_Plot[N_∑T_Plot], ymin=0.0, xtick=(DateTick,DateTick2), xrotation=rad2deg(pi/4), framestyle=:box, grid=true)

				end # if: option.hyPix.Plot_θ
				
				Plots.savefig(Plot, Path)
				println("			 ~ ", Path, "~")
			
				return nothing
				end  # function: TIMESERIES

			end  # module: plots
			# ............................................................

end  # module plotHypix
# ............................................................