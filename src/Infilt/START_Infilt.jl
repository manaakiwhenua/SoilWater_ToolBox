# =============================================================
#		MODULE: infiltration
# =============================================================
module infilt
	import ..option, ..sorptivity, ..param, ..wrc, ..kunsat, ..option, ..infiltInitialize, ..bestFunc, ..stats, ..tool, ..quasiExact, ..ofBest, ..hydroRelation
	import BlackBoxOptim, Statistics
	export START_INFILTRATION

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : START_INFILT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function START_INFILTRATION(∑Infilt_Obs, ∑Psd, hydro, hydroInfilt, Id_Select, infiltParam, N_Infilt, N_SoilSelect, Tinfilt)

		# INITIALIZE
			T, infiltOutput, hydroInfilt, ∑Infilt_3D, ∑Infilt_1D = infiltInitialize.INFILT_INITIALIZE(∑Infilt_Obs, ∑Psd, hydroInfilt, infiltParam, N_Infilt, N_SoilSelect, Tinfilt)

		for iZ=1:N_SoilSelect
			println( "iZ= $iZ")

			# No optimization required running from hydro derived from laboratory
			if option.infilt.OptimizeRun == :Run && option.globalopt.θΨ ≠ :No #<>=<>=<>=<>=<>
				# Hydraulic param from laboratory
					hydroInfilt = deepcopy(hydro)
				 
				# Not to have errors
					infiltParam.θ_Ini[iZ] = max(hydroInfilt.θr[iZ] + eps(), infiltParam.θ_Ini[iZ])

				if option.infilt.Model == :Best_Univ
					∑Infilt_3D, T_TransStead =  bestFunc.BEST_UNIVERSAL_START(∑Infilt_3D, hydroInfilt, infiltOutput, infiltParam, iZ, N_Infilt, T)

				elseif option.infilt.Model == :QuasiExact
					∑Infilt_3D = quasiExact.HYDRO_2_INFILTRATION3D(∑Infilt_3D, hydroInfilt, infiltParam, iZ, N_Infilt, T)

				end # option.infilt.Model


			elseif option.infilt.OptimizeRun == :RunOptKs && option.globalopt.θΨ ≠ :No #<>=<>=<>=<>=<>	
				# Hydraulic param from laboratory
					hydroInfilt = deepcopy(hydro)
				
				# Not to have errors
				infiltParam.θ_Ini[iZ] = max(hydroInfilt.θr[iZ] + eps(), infiltParam.θ_Ini[iZ])
				
				SearchRange =[(log10(hydroInfilt.Ks_Min), log10(hydroInfilt.Ks_Max))]

				Optimization = BlackBoxOptim.bboptimize(P -> OF_INFILT_2_HYDRO(∑Infilt_3D, ∑Infilt_Obs, hydroInfilt, infiltOutput, infiltParam, iZ, N_Infilt, T; Ks=10.0^P[1])[1]; SearchRange=SearchRange, NumDimensions=1, TraceMode=:silent)

				hydroInfilt.Ks[iZ] = 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[1]


			elseif option.infilt.OptimizeRun == :Opt && option.infilt.HydroModel == :Kosugi # <>=<>=<>=<>=<>	

				if option.infilt.σ_2_Ψm == :Constrained
					SearchRange =[ (hydroInfilt.σ_Min[iZ], hydroInfilt.σ_Max[iZ]), (0.0, 1.0), (log10(hydroInfilt.Ks_Min[iZ]), log10(hydroInfilt.Ks_Max[iZ]))]

					Optimization = BlackBoxOptim.bboptimize(P -> OF_INFILT_2_HYDRO(∑Infilt_3D, ∑Infilt_Obs, hydroInfilt, infiltOutput, infiltParam, iZ, N_Infilt, T; σ=P[1], Ψm=P[2], Ks=10.0^P[3])[1]; SearchRange=SearchRange, NumDimensions=3, TraceMode=:silent)

					hydroInfilt.σ[iZ]  = BlackBoxOptim.best_candidate(Optimization)[1]
					hydroInfilt.Ψm[iZ] = BlackBoxOptim.best_candidate(Optimization)[2]
					hydroInfilt.Ks[iZ] = 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[3]
		
					hydroInfilt = hydroRelation.FUNCTION_σ_2_Ψm_SOFTWARE(hydroInfilt, iZ, option.infilt)
						
				else
					SearchRange =[ (hydroInfilt.σ_Min[iZ], hydroInfilt.σ_Max[iZ]), (log10(hydroInfilt.Ψm_Min[iZ]), log10(hydroInfilt.Ψm_Max[iZ])), (log10(hydroInfilt.Ks_Min[iZ]), log10(hydroInfilt.Ks_Max[iZ]))]

					Optimization = BlackBoxOptim.bboptimize(P -> OF_INFILT_2_HYDRO(∑Infilt_3D, ∑Infilt_Obs, hydroInfilt, infiltOutput, infiltParam, iZ, N_Infilt, T; σ=P[1], Ψm=10.0^P[2], Ks=10.0^P[3])[1]; SearchRange=SearchRange, NumDimensions=3, TraceMode=:silent)

					hydroInfilt.σ[iZ]  = BlackBoxOptim.best_candidate(Optimization)[1]
					hydroInfilt.Ψm[iZ] = 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[2]
					hydroInfilt.Ks[iZ] = 10.0 ^ BlackBoxOptim.best_candidate(Optimization)[3]
				end # option.infilt.σ_2_Ψm
	
			else
				error("ERROR SoilWaterToolBox = $(option.infilt.Model) not found")
			end # option.infilt.OptimizeRun


			# OUTPUTS RUNNING THE OPTIMAL INFILTRATION
				infiltOutput.Sorptivity[iZ] = sorptivity.SORPTIVITY(infiltParam.θ_Ini[iZ], iZ, hydroInfilt) 

				if option.infilt.Model == :Best_Univ 
					∑Infilt_3D, T_TransStead = bestFunc.BEST_UNIVERSAL_START(∑Infilt_3D, hydroInfilt, infiltOutput, infiltParam, iZ, N_Infilt, T)

				elseif option.infilt.Model == :QuasiExact
					# ∑Infilt_3D = quasiExact.HYDRO_2_INFILTRATION3D(∑Infilt_3D, hydroInfilt, infiltParam, iZ, N_Infilt, T)
					∑Infilt_3D = quasiExact.HYDRO_2_INFILTRATION3D(∑Infilt_3D, hydroInfilt, infiltParam, iZ, N_Infilt, T)

				end # option.infilt.Model

			# Statistics
				iT_TransSteady = infiltOutput.iT_TransSteady_Data[iZ]

				infiltOutput.Nse_Trans[iZ] = 1.0 - stats.NASH_SUTCLIFE_MINIMIZE( ∑Infilt_Obs[iZ,2:iT_TransSteady], ∑Infilt_3D[iZ, 2:iT_TransSteady]; Power=2.0)

				infiltOutput.Nse_Steady[iZ] = 1.0 - stats.NASH_SUTCLIFE_MINIMIZE( log10.(∑Infilt_Obs[iZ,iT_TransSteady+1:N_Infilt[iZ]]), log10.(∑Infilt_3D[iZ,iT_TransSteady+1:N_Infilt[iZ]]); Power=2.0)

				infiltOutput.Nse[iZ] = 0.5 * infiltOutput.Nse_Trans[iZ] + 0.5 * infiltOutput.Nse_Steady[iZ]

		end # for iZ=1:N_SoilSelect


		# AVERAGE STATISTICS
         Nse_Trans  = Statistics.mean(infiltOutput.Nse_Trans[1:N_SoilSelect])
         Nse_Steady = Statistics.mean(infiltOutput.Nse_Steady[1:N_SoilSelect])
         Nse        = Statistics.mean(infiltOutput.Nse[1:N_SoilSelect])

			println("    ~ $(option.infilt.SorptivityModel) ~")
			println("    ~ Nse= $Nse Nse_Trans= $Nse_Trans,  Nse_Steady= $Nse_Steady ~")


		# CONVERTING INFILTRATION DIMENSIONS
			if option.infilt.Model == :Best_Univ
				for iZ=1:N_SoilSelect
					∑Infilt_1D = bestFunc.CONVERT_3D_2_1D(∑Infilt_3D, ∑Infilt_1D, hydroInfilt, infiltParam, iZ, N_Infilt, T)
				end
			elseif option.infilt.Model == :QuasiExact
				for iZ=1:N_SoilSelect
					∑Infilt_1D = quasiExact.CONVERT_3D_2_1D(∑Infilt_3D, ∑Infilt_1D, hydroInfilt, infiltParam, iZ, N_Infilt, T)
				end	
			end #  option.infilt

		return infiltOutput, hydroInfilt, ∑Infilt_3D, ∑Infilt_1D

	end # FUNCTION: START_INFILTRATION


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OF_INFILT_2_HYDRO
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function OF_INFILT_2_HYDRO(∑Infilt_3D, ∑Infilt_Obs, hydroInfilt, infiltOutput, infiltParam, iZ, N_Infilt, T; σ=hydroInfilt.σ[iZ], Ψm=hydroInfilt.Ψm[iZ], θr=hydroInfilt.θr[iZ], θs=hydroInfilt.θs[iZ], Ks=hydroInfilt.Ks[iZ])

         hydroInfilt.θs[iZ]       = θs
         hydroInfilt.θr[iZ]       = θr
         hydroInfilt.Ks[iZ]       = Ks
         hydroInfilt.σ[iZ]        = σ
         hydroInfilt.Ψm[iZ]       = Ψm
         hydroInfilt.θsMacMat[iZ] = θs
         hydroInfilt.ΨmMac[iZ]    = Ψm
         hydroInfilt.σMac[iZ]     = σ

			hydroInfilt = hydroRelation.FUNCTION_σ_2_Ψm_SOFTWARE(hydroInfilt, iZ, option.infilt)

			if option.infilt.Model == :Best_Univ
				return Nse = ofBest.OF_BEST(∑Infilt_3D, ∑Infilt_Obs, hydroInfilt, infiltOutput, infiltParam, iZ, N_Infilt, T; W=0.5)

			elseif option.infilt.Model == :QuasiExact
				quasiExact.OF_QUASIEXACT(∑Infilt_Obs, hydroInfilt, infiltOutput, infiltParam, iZ, N_Infilt, T)
			end # Option.infilt.Model

		end # FUNCTION: OF_INFILT_2_HYDRO_Best


end # MODULE: infilt
	