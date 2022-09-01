# =============================================================
#		MODULE: hydro
# =============================================================
module hypix

	import ..climate, ..discretization, ..horizonLayer, ..hydroStruct, ..hypixModel, ..hypixOpt, ..interpolate, ..memory, ..ofHypix, ..option, ..param, ..pathHypix, ..plotHypix, ..plotOther, ..readHypix, ..reading, ..stats, ..tableHypix, ..thetaObs, ..tool, ..vegStruct, ..waterBalance, ..Δtchange, ..θaver, ..sitename
	
	import Statistics, Dates, DelimitedFiles
	export HYPIX_START

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYPIX_START
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function HYPIX_START()

		# ===========================================================
		# 					LOOP FOR DIFFERENTY SIMULATIONS
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		for iSim = param.hyPix.iSim_Start : param.hyPix.iSim_End
			println("=== START RUNNING Hypix_1D ==== \n")
			println("		~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
			println("		=== === === START: Looping with time $iSim steps \n")


		# READING STRUCTURE OF PATH
			pathHyPix = pathHypix.PATHHYPIX(iSim)

				println("\n			==== ==== ===  $(pathHyPix.SiteName_Hypix) 	=== ==== ====\n")

			# COUNT SIMULATIONS
				iSim_Count = iSim - param.hyPix.iSim_Start + 1

			# ACTUAL TIME
				Time_Start = Dates.now()

			# DISCRETISATION ~~~~~
				# Read discretisaation
					Layer, N_iHorizon, N_iZ, Z, θ_Ini = readHypix.DISCRETIZATION(pathHyPix)
				# Process discretisation ~~~~~
					discret = discretization.DISCRETIZATION(N_iZ, Z)

			# CLIMATE DATA  ~~~~~
				# Read climate
					clim = readHypix.CLIMATE(pathHyPix)
				# Process climate
					∑Pet_Climate, ∑Pr_Climate, ∑T_Climate, N_∑T_Climate, Temp = climate.CLIMATE(clim)

			# OBSERVED θ  ~~~~~
			if option.hyPix.θobs
				# Read observed θ
					obsθ = readHypix.TIME_SERIES(pathHyPix)
				
				# Process observed θ
					obsθ = thetaObs.ΘOBS(obsθ, clim, discret, Z)
			end #  option.hyPix.θobs 

			# MEMORY
				if iSim_Count == 1
					N_∑T_Plot = param.hyPix.iSim_End - param.hyPix.iSim_Start + 1

					global Efficiency                 = fill(0.0::Float64, N_∑T_Plot)
					global Global_WaterBalance        = fill(0.0::Float64, N_∑T_Plot)
					global Global_WaterBalance_NormPr = fill(0.0::Float64, N_∑T_Plot)
					global RmseBest                   = fill(0.0::Float64, N_∑T_Plot)
					global SwcRoots                   = fill(0.0::Float64, N_∑T_Plot)
					global WofBest                    = fill(0.0::Float64, N_∑T_Plot)
					global ΔRunTimeHypix              = fill(0.0::Float64, N_∑T_Plot)
					global ΔT_Average                 = fill(0.0::Float64, N_∑T_Plot)
					global ∑ΔQ_Bot                    = fill(0.0::Float64, N_∑T_Plot)
					global ∑∑ΔSink                    = fill(0.0::Float64, N_∑T_Plot)
				end
				∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑Pet, ∑Pr, ∑T, CropCoeficientᵀ, iNonConverge_iSim, Laiᵀ, N_Memory, Q, Residual, ΔEvaporation, ΔHpond, ΔΨmax, ΔPet, ΔPr, ΔSink, ΔT, θ, θSim, Ψ, Ψ_Max, Ψ_Min, Ψbest = memory.MEMORY(clim, iSim_Count, N_∑T_Climate, N_iZ, obsθ)

			# INITIALIZING THE STRUCTURE
				# Initializing hydroHorizon structure
					hydroHorizon = hydroStruct.HYDROSTRUCT(N_iHorizon)

				# Initializing hydraulic param into structure 
					hydro = hydroStruct.HYDROSTRUCT(N_iZ)

				# Initialiozing  vegetation parameters into veg structure
					veg = vegStruct.VEGSTRUCT()

				# Optimisation
					hydroHorizon_best = hydroStruct.HYDROSTRUCT(N_iHorizon)
					hydro_best        = hydroStruct.HYDROSTRUCT(N_iZ)
					veg_best          = vegStruct.VEGSTRUCT()
				
			# ========================================================
	 
			# OBTAINING HYDRAULIC AND VEGETATION PARAMETERS
			if option.hyPix.Optimisation
				hydro, hydroHorizon, optim, veg = readHypix.HYPIX_PARAM(Layer, hydro, hydroHorizon, iSim, N_iZ, veg)
			else
			# Reading veg parameters
				veg, ~ = reading.READ_STRUCT(veg, pathHyPix.HyPix_VegParam)
	
			# Hydraulic parameters
				hydroHorizon, ~ = reading.READ_STRUCT(hydroHorizon, pathHyPix.HyPix_HydroParam)

				hydro        = horizonLayer.HYDROHORIZON_2_HYDRO(N_iZ, Layer, hydroHorizon)

			# options of optim
			
				Flag_Opt = false
				NparamOpt = 0

				optim = ( NparamOpt=NparamOpt, Flag_Opt=Flag_Opt)		
			end # option.hyPix.Optimisation

			# SINK TERM 
            Laiᵀ_η            = readHypix.LOOKUPTABLE_LAI(clim, pathHyPix, veg)
            CropCoeficientᵀ_η = readHypix.LOOKUPTABLE_CROPCOEFICIENT(clim, pathHyPix, veg)

			if optim.Flag_Opt
				hydro, hydro_best, hydroHorizon, hydroHorizon_best, veg, veg_best, WofBest = hypixOpt.HYPIXOPTIMISATION_START(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑Pet, ∑Pet_Climate, ∑Pr, ∑Pr_Climate, ∑T, ∑T_Climate, obsθ, clim, CropCoeficientᵀ, CropCoeficientᵀ_η, discret, Layer, hydro, hydro_best, hydroHorizon, hydroHorizon_best, iSim_Count, Laiᵀ, Laiᵀ_η, N_∑T_Climate, N_iHorizon, N_iZ, optim, Q, Residual, veg, veg_best, WofBest, Z, ΔEvaporation, ΔHpond, ΔΨmax, ΔPet, ΔPr, ΔSink, ΔT, θ, θ_Ini, θSim, Ψ, Ψ_Max, Ψ_Min, Ψbest)
			end
			
			# if Flag_Opt then it will rerun with the optimal parameters
			∑Pet, ∑Pr, ∑T, ∑T_Climate, clim, discret, iNonConverge, Iter_CountTotal, N_iRoot, N_iT, N_iZ, Q, veg, ΔEvaporation, ΔHpond, ΔRootDensity, ΔT, θ, Ψ = hypixModel.HYPIX(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, ∑Pet, ∑Pet_Climate, ∑Pr, ∑Pr_Climate, ∑T, ∑T_Climate, clim, CropCoeficientᵀ, CropCoeficientᵀ_η, discret, hydro, Laiᵀ, Laiᵀ_η, N_∑T_Climate, N_iZ, Q, Residual, veg, Z, ΔEvaporation, ΔHpond, ΔΨmax, ΔPet, ΔPr, ΔSink, ΔT, θ, θ_Ini, Ψ, Ψ_Max, Ψ_Min, Ψbest)

			# WATER BALANCE
				# Computed after the warmup period
				∑∑WaterBalance, ∑WaterBalance_η, ∑ΔSink, i∑T_CalibrStart, ΔStorage = waterBalance.WATERBALANCE(∑T, obsθ, discret, hydro, N_iRoot, N_iT, N_iZ, Q, ΔSink, ΔT, θ, Ψ)

			# SIGNATURE
			# if option.hyPix.Signature_Run
			# 	Signature_Deficit_Obs, Signature_Max_Obs, Signature_Saturated_Obs, Signature_Senescence_Obs, Signature_Deficit_Sim, Signature_Max_Sim, Signature_Saturated_Sim, Signature_Senescence_Sim = signature.SIGNATURE(∑T[1:N_iT], obsθ, hydroHorizon, N_iRoot, N_iT, N_iZ, veg, ΔRootDensity, Ψ[1:N_iT,1:N_iZ])
			# end 

			# SUMMARY HOW GOOD THE SIMULATION
				# Computed climate day after the warmup period
				i∑T_CalibrStart_Day = 1::Int64 
				while ∑T_Climate[i∑T_CalibrStart_Day] < obsθ.∑T[1]
					i∑T_CalibrStart_Day += 1
				end
				i∑T_CalibrStart_Day += 1
				∑Pr_Clim = sum(clim.Pr[i∑T_CalibrStart_Day:clim.N_Climate]) 
				∑Pet_Net = sum(clim.Pet[i∑T_CalibrStart_Day:clim.N_Climate])

				# Soil water content for the rootzone at the end of simulation
				@fastmath @inbounds @simd for iZ=1:N_iRoot
					SwcRoots[iSim_Count] += θ[N_iT,iZ] * discret.ΔZ[iZ]
				end

				# Timing 
					Time_End = Dates.now()
					ΔRunTimeHypix[iSim_Count] = Dates.value(Time_End - Time_Start) / 1000

				# Convergence rate
					iNonConverge_iSim[iSim_Count]          = iNonConverge

					Efficiency[iSim_Count]                 = ceil(Int, 86400 * Iter_CountTotal / ∑T[N_iT])
					
					ΔT_Average[iSim_Count]                 = ceil(Int, Statistics.mean(ΔT[i∑T_CalibrStart:N_iT]))
					
					Global_WaterBalance[iSim_Count]        = ∑∑WaterBalance
					
					Global_WaterBalance_NormPr[iSim_Count] = 100.0 * ∑WaterBalance_η[N_iT]
					
					∑∑ΔSink[iSim_Count]                    = ∑ΔSink[N_iT]
				
				# Ground water recharge
					∑ΔQ_Bot[iSim_Count] = 0.0
					@fastmath @inbounds @simd for iT=i∑T_CalibrStart:N_iT
						∑ΔQ_Bot[iSim_Count] += ΔT[iT] * Q[iT, N_iZ+1]
					end
           
				println("		=== ===START: summary  $iSim steps ...")

					println("			∑Pr 			= ", ceil(Int, ∑Pr_Clim), "  [mm]")
					println("			∑Pr_Soil 		= ", ceil(Int, ∑Pr[N_iT] - ∑Pr[i∑T_CalibrStart]),  "  [mm]")
					println("			∑Pr_Intercepted/∑Pr 	= ", ceil(Int, 100. * (∑Pr_Clim - (∑Pr[N_iT]-∑Pr[i∑T_CalibrStart])) / ∑Pr_Clim),  "  [%]")
					println("			∑Pet_Net 		= ", ceil(Int, ∑Pet_Net), "  [mm]")
					println("			∑Pet 			= ", ceil(Int, ∑Pet[N_iT]- ∑Pet[i∑T_CalibrStart]), "  [mm]")
					println("			∑ΔSink/∑Pet_Net 	= ", ceil(Int, 100.0 * ∑ΔSink[N_iT] / ∑Pet_Net), "  [%] \n")
					
					println("			∑SoilWaterContentRootEnd = ", ceil(Int, SwcRoots[iSim_Count]), "  [mm]")
					println("			∑ΔSink 			= ", -ceil(Int, ∑∑ΔSink[iSim_Count]), "  [mm]")
					println("			∑Infilt_Bot 		= ", -ceil( Int, ∑ΔQ_Bot[iSim_Count] ), "  [mm]")
					println("			ΔHpond at end 		= ", ceil(Int, ΔHpond[N_iT]), "  [mm] \n")

					println("			iNonConverge 			= ", iNonConverge_iSim[iSim_Count], "  [count]")
					println("			Global_WaterBalance_NormPr 	= ", round(Global_WaterBalance_NormPr[iSim_Count], digits=2), "  [%]")
					println("			Global_WaterBalance 		= ", 	round(Global_WaterBalance[iSim_Count], digits=2), "  [mm]")
					println("			Average ΔT 			= ",  ΔT_Average[iSim_Count] , "  [seconds]")
					println("			ΔTmin 				= ",   round(minimum(ΔT[i∑T_CalibrStart:N_iT]), digits=0) , "  [seconds]")
					println("			ΔTmax 				= ",  round(maximum(ΔT[i∑T_CalibrStart:N_iT]), digits=0) , "  [seconds]")
					println("			ΔT_HyPix 			= ", ceil(Int, ΔRunTimeHypix[iSim_Count]) , "  [seconds]")			
					println("			Efficiency 			= ", Efficiency[iSim_Count], "  [iTer day-1], \n")


				∑T_Plot, ∑T_Date_Plot, ∑WaterBalance_η_Plot, Date_Plot, Flag_Plot_Pond, N_∑T_Plot, ΔEvaporation_Plot, ΔFlux_Plot, ΔPet_Plot, ΔPond_Plot, ΔPr_Plot, ΔSink_Plot, ΔT_Plot, θ_Plot, θobs_Plot, Ψ_Plot = Δtchange.CHANGE_OUTPUT_ΔT(∑Pet[1:N_iT], ∑Pr[1:N_iT], ∑T[1:N_iT], ∑WaterBalance_η[1:N_iT], ∑ΔSink[1:N_iT], obsθ, clim, N_iT, N_iZ, Q[1:N_iT,1:N_iZ+1], veg, ΔEvaporation[1:N_iT], ΔHpond[1:N_iT], ΔT[1:N_iT], θ[1:N_iT,1:N_iZ], Ψ[1:N_iT,1:N_iZ], ∑T_Climate, pathHyPix)

			# Computing average simulated θ to comapre it with average observed θ
			if option.hyPix.θobs_Average && option.hyPix.θobs	
				θsim_Aver = θaver.θAVER(discret; Z=Z, θ_Plot=θ_Plot, N_iZ=N_iZ, N_∑T_Plot=N_∑T_Plot, Zaver=400.0)

				RmseBest[iSim_Count] = stats.NSE(θobs_Plot[1:N_∑T_Plot], θsim_Aver[1:N_∑T_Plot])
				
			elseif !(option.hyPix.θobs) && option.hyPix.θobs
				RmseBest[iSim_Count] = ofHypix.θof.RMSE_θ(∑T, obsθ, N_iT, N_iZ, θ, θSim)
				θsim_Aver = []
			end

			if option.hyPix.θobs
				println("			Nse 			= ", round(RmseBest[iSim_Count], digits=5), "  [mm3 mm-3]")
			end
			println("		=== === END: summary \n")


			if option.hyPix.Table
			println("		=== === START: Table === ===")
	
				# Writing values of hydraulic parameters
				tableHypix.HYDRO(hydroHorizon, iSim, N_iHorizon, pathHyPix)

				# Writing values of veg parameters
				tableHypix.VEG(veg, iSim, pathHyPix)

				SiteNames = sitename.SITENEAME()

				tableHypix.PERFORMANCE(∑∑ΔSink, ∑ΔQ_Bot, Efficiency, Global_WaterBalance, Global_WaterBalance_NormPr, iNonConverge_iSim, iSim, RmseBest, SwcRoots, WofBest, ΔRunTimeHypix, ΔT_Average, SiteNames[param.hyPix.iSim_Start:param.hyPix.iSim_End], pathHyPix)		

				if option.hyPix.Table_Discretization
					tableHypix.DISCRETIZATION(discret, N_iZ, Z[1:N_iZ], pathHyPix)
				end
				if  option.hyPix.Table_TimeSeries
					tableHypix.TIME_SERIES(∑T[1:N_iT], ΔT[1:N_iT], ∑Pr[1:N_iT], ΔPr[1:N_iT], ΔHpond[1:N_iT], ΔT[1:N_iT].*Q[1:N_iT,1], ∑WaterBalance_η[1:N_iT], iSim, pathHyPix)
				end
				if option.hyPix.Table_TimeSeriesDaily
					tableHypix.TIME_SERIES_DAILY(∑T_Plot[1:N_∑T_Plot], ∑WaterBalance_η_Plot[1:N_∑T_Plot], Date_Plot[1:N_∑T_Plot], iSim, N_∑T_Plot, ΔEvaporation_Plot[1:N_∑T_Plot], ΔFlux_Plot[1:N_∑T_Plot, N_iZ+1], ΔPet_Plot[1:N_∑T_Plot], ΔPond_Plot[1:N_∑T_Plot], ΔPr_Plot[1:N_∑T_Plot], ΔSink_Plot[1:N_∑T_Plot], pathHyPix)
				end
				if option.hyPix.Table_θ
					tableHypix.θ(∑T[1:N_iT], θ[1:N_iT,1:N_iZ], discret.Znode[1:N_iZ], iSim, pathHyPix)
				end
				if option.hyPix.Table_Ψ
					tableHypix.Ψ(∑T[1:N_iT], Ψ[1:N_iT,1:N_iZ], discret.Znode[1:N_iZ], iSim, pathHyPix)
				end
				if option.hyPix.Table_Q
					tableHypix.Q(∑T[1:N_iT], Q[1:N_iT,1:N_iZ+1], Z[N_iZ], discret.Znode[1:N_iZ], iSim, pathHyPix)
				end
				if option.hyPix.Tabule_θΨ
					tableHypix.θΨ(hydroHorizon, iSim, N_iHorizon, pathHyPix)
					tableHypix.KΨ(hydroHorizon, iSim, N_iHorizon, pathHyPix)
				end
				if option.hyPix.Table_Climate
					tableHypix.DAILY_CLIMATE(∑T_Climate, clim, iSim, pathHyPix)
				end
				if option.hyPix.θobs_Average && option.hyPix.θobs
					tableHypix.θAVERAGE(Date_Plot[1:N_∑T_Plot], iSim, θobs_Plot[1:N_∑T_Plot], θsim_Aver[1:N_∑T_Plot], pathHyPix)
				end
				# if option.hyPix.Signature_Run
				# 	tableHypix.SIGNATURE(iSim, Signature_Deficit_Obs, Signature_Max_Obs, Signature_Saturated_Obs, Signature_Senescence_Obs, Signature_Deficit_Sim, Signature_Max_Sim, Signature_Saturated_Sim, Signature_Senescence_Sim)
				# end
			println("		=== === END: Table === === \n")
			end  # if option.hyPix.Table
				
			if option.globalopt.Ploting
			println("		=== === START: Plotting === ===")

				if option.hyPix.Plot_Other
					# plotOther.ΨMINΨMAX(hydro, pathHyPix)
					plotOther.WOF_STEPS(pathHyPix)
					# plotOther.SE_Ψ_CONSTRAINED(hydro, pathHyPix)
					# plotOther.PLOT_σ_2_θr(hydro, pathHyPix)
					# plotOther.PLOT_θΨ_Δθ(hydro, pathHyPix)
				end # option.hyPix.Plot_Other
				
				if option.hyPix.Plot_Hypix
					plotHypix.plots.TIMESERIES(∑T_Date_Plot, ∑T_Plot, obsθ, discret, Flag_Plot_Pond, iSim, N_∑T_Plot, N_iZ, ΔEvaporation_Plot, ΔFlux_Plot, ΔPet_Plot, ΔPond_Plot, ΔPr_Plot, ΔSink_Plot, θ_Plot, θobs_Plot, clim, i∑T_CalibrStart_Day, θsim_Aver, pathHyPix)
					# plotHypix.TIME_SERIES(∑T_Plot, ∑WaterBalance_η_Plot, obsθ, discret, Flag_Plot_Pond, iSim, N_∑T_Plot, N_iZ, ΔFlux_Plot, ΔPet_Plot, ΔPond_Plot, ΔPr_Plot, ΔSink_Plot, ΔT_Plot, θ_Plot, θobs_Plot, Ψ_Plot, pathHyPix)
				end
				if option.hyPix.Plot_θΨK
					plotHypix.θΨK(hydroHorizon, N_iHorizon, iSim, pathHyPix)
				end
				if option.hyPix.Plot_Vegetation && option.hyPix.RootWaterUptake
					plotHypix.VEG_FUNCTIONS(discret, iSim, N_iRoot, veg, Z, ΔRootDensity, pathHyPix)
				end
				if option.hyPix.Plot_Interception
					plotHypix.plots.RAINFALL_INTERCEPTION(clim, i∑T_CalibrStart_Day, iSim, pathHyPix)
				end
				if  option.hyPix.Plot_Sorptivity
					plotHypix.plots.PLOT_SORPTIVITY(iSim, hydro, pathHyPix)
				end
			println("		=== === END: Plotting === === \n")
			end # if option.hyPix.Plotting
	
	println("	=== === === END  : Looping with time ")
	println("	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n \n")

	# Garbage collect
	GC.gc()
	end # for loop: iSim

	println("\n ==== END RUNNING Hypix_1D =====")

	return
	end # function HYPIX_START
	
end  # module hydro
# ............................................................