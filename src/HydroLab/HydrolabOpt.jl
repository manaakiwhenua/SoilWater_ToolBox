# =============================================================
#		module: hypixOpt
# =============================================================
module hydrolabOpt

	import ..OfHydrolab, ..option, ..param, ..tool, ..optimize, ..hydroRelation, ..psdThetar
	using BlackBoxOptim, Statistics
	export HYDROLABOPT_START

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYPIXOPT_START
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function HYDROLABOPT_START(;N_SoilSelect, ∑Psd, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ=[0], Ψ_KΨ=[0], N_KΨ=1, hydro, hydroOther, optionHydro, optim)

		for iZ = 1:N_SoilSelect
			# CORRECTION OF THE FEASIBLE RANGE ~~~
				θ_Min = minimum(θ_θΨ[iZ, 1:N_θΨ[iZ]])  	# Smallest measure θ

				θ_Max = maximum(θ_θΨ[iZ, 1:N_θΨ[iZ]])  	# Greatest measure θ

			# CORRECTING Θr ~~~
				if ("θr" ∈ optim.ParamOpt)
					hydro.θr_Max[iZ] = max( min(θ_Min-0.005, hydro.θr_Max[iZ]), hydro.θr_Min[iZ]) # Maximum value of θr

					# Changing the feasible range of θr
					iθr = findfirst(isequal("θr"), optim.ParamOpt)[1]
					optim.ParamOpt_Max[iθr] = hydro.θr_Max[iZ]

				elseif ("θr" ∉ optim.ParamOpt) && (optionHydro.θrOpt==:ParamPsd) && (option.globalopt.Psd) # Derive θr frpm PSD
					hydro.θr[iZ] = min(psdThetar.PSD_2_θr_FUNC(∑Psd, hydro, iZ), θ_Min-0.005)

				end # if ("θr" ∈ optim.ParamOpt)

			# TEST IF EXIST Ψ=0
				if minimum(Ψ_θΨ[iZ,1:N_θΨ[iZ]]) < eps(1000.0)
					Flag_Ψ0 = true
				else
					Flag_Ψ0 = false
				end

			# CORRECTING θS  ~~~
				if ("θs" ∈ optim.ParamOpt) && Flag_Ψ0
					hydro.θs_Min[iZ] = θ_Max * 0.75
					hydro.θs_Max[iZ] = θ_Max
					hydro.Φ[iZ] = θ_Max / 0.95

					# Modifying the searchrange
						iθs = findfirst(isequal("θs"), optim.ParamOpt)[1]
						optim.ParamOpt_Min[iθs] = hydro.θs_Min[iZ]
						optim.ParamOpt_Max[iθs] = hydro.θs_Max[iZ]

				elseif ("θs" ∉ optim.ParamOpt) && Flag_Ψ0 # <>=<>=<>=<>=<>
						hydro.θs[iZ] = θ_Max
						hydro.Φ[iZ] = hydro.θs[iZ] / 0.95

				elseif  ("θs" ∉ optim.ParamOpt) && (optionHydro.θsOpt == :Φ) # <>=<>=<>=<>=<>
						if hydro.Φ[iZ] * 0.95 > θ_Max + 0.005
							hydro.θs[iZ] = hydro.Φ[iZ] * 0.95
						elseif hydro.Φ[iZ] * 0.965 > θ_Max + 0.005
							hydro.θs[iZ] = hydro.Φ[iZ] * 0.965
						else
							hydro.θs[iZ] = max(hydro.Φ[iZ] - 0.005, hydro.θs_Max[iZ] + 0.005)
						end # hydro.Φ[iZ] * 0.95 > hydro.θs_Min[iZ]

				elseif ("θs" ∈ optim.ParamOpt) && (option.hydro.θs_MinFromData )# <>=<>=<>=<>=<>
					hydro.θs_Min[iZ] = θ_Max + 0.005
					hydro.θs_Max[iZ] = max(θ_Max + 0.25, hydro.θs_Max[iZ])

					# Modifying the searchrange
						iθs = findfirst(isequal("θs"), optim.ParamOpt)[1]
						optim.ParamOpt_Min[iθs] = hydro.θs_Min[iZ]
						optim.ParamOpt_Max[iθs] = hydro.θs_Max[iZ]
		
				end # optionHydro.θsOpt
				
			# CORRECTING Ks

				# TEST IF EXIST Ψ=0
				if ("Ks" ∈ optim.ParamOpt)
					if minimum(Ψ_KΨ[iZ,1:N_θΨ[iZ]]) < eps(1000.0)
						Flag_K0 = true
					else
						Flag_K0 = false
					end
				end

				if ("Ks" ∈ optim.ParamOpt) && Flag_K0
					hydro.Ks_Max[iZ] = maximum(K_KΨ[iZ, 1:N_KΨ[iZ]]) # Greatest measure of Kunsat)

					# Modifying the searchrange
						iKs = findfirst(isequal("Ks"), optim.ParamOpt)[1]
						optim.ParamOpt_Max[iKs] = hydro.Ks_Max[iZ]

				elseif ("Ks" ∈ optim.ParamOpt) && option.hydro.Ks_MinMaxFromData
					Ks_max = maximum(K_KΨ[iZ, 1:N_KΨ[iZ]])
					hydro.Ks_Min[iZ] = Ks_max # Greatest measure of Kunsat
					hydro.Ks_Max[iZ] = max(hydro.Ks_Min[iZ] + 0.01, hydro.Ks_Max[iZ])

					# Modifying the searchrange
						iKs = findfirst(isequal("Ks"), optim.ParamOpt)[1]
						optim.ParamOpt_Min[iKs] = hydro.Ks_Min[iZ]
						optim.ParamOpt_Max[iKs] = hydro.Ks_Max[iZ]

				elseif ("Ks" ∉ optim.ParamOpt) && (optionHydro.KsOpt == :Data) && (option.hydro.KunsatΨ)
					Ks_max = maximum(K_KΨ[iZ, 1:N_KΨ[iZ]])
					hydro.Ks[iZ] = Ks_max

				end # "Ks" ∈ optim.ParamOpt
			
			# Updated searchrange
				SearchRange = optimize.SEARCHRANGE(optim)

			# OPTIMIZATION: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

				Optimization = BlackBoxOptim.bboptimize(X -> hydrolabOpt.OF_HYDROLAB(X, hydro, iZ, K_KΨ, N_KΨ, N_θΨ, optim, optionHydro, θ_θΨ, Ψ_KΨ, Ψ_θΨ); SearchRange=SearchRange, NumDimensions=optim.NparamOpt, TraceMode=:silent)

				X = BlackBoxOptim.best_candidate(Optimization)

				hydro = hydrolabOpt.PARAM_2_hydro(hydro, iZ, optim, optionHydro, X)

				# STATISTICS
					Of, Of_θΨ, Of_Kunsat = OfHydrolab.OF_WRC_KUNSAT(iZ, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro, optim, optionHydro) 

					hydroOther.Rmse[iZ], hydroOther.Rmse_KΨ[iZ], hydroOther.Rmse_θΨ[iZ] = OfHydrolab.OF_RMSE(iZ, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro, optim, optionHydro) 

					hydroOther.Nse[iZ]    = 1.0 - Of
					hydroOther.Nse_θΨ[iZ] = 1.0 - Of_θΨ

					if optionHydro.KunsatΨ
						hydroOther.Nse_KΨ[iZ] = 1.0 - Of_Kunsat
					end
		end # for iZ = 1:N_SoilSelect

		# OVERALL STATISTICS OF THE OPTIMIZATION
			Nse_θΨ_Aver = Statistics.mean(hydroOther.Nse_θΨ[1:N_SoilSelect])
			Nse_KΨ_Aver = Statistics.mean(max.(hydroOther.Nse_KΨ[1:N_SoilSelect],0.0))

			Rmse_Aver    = Statistics.mean(hydroOther.Rmse[1:N_SoilSelect])
			Rmse_θΨ_Aver = Statistics.mean(hydroOther.Rmse_θΨ[1:N_SoilSelect])
			Rmse_KΨ_Aver = Statistics.mean(hydroOther.Rmse_KΨ[1:N_SoilSelect])
				
			if option.hydro.KunsatΨ
				Nse_Aver = (Nse_θΨ_Aver + Nse_KΨ_Aver) / 2.0
			else
				Nse_Aver = Nse_θΨ_Aver
			end
			
			println("\n    ==  Optimizing Hydraulic parameters == ")
			println("    	~  Nse_θΨ = $(round(Nse_θΨ_Aver,digits=3)),  Nse_KΨ = $(round(Nse_KΨ_Aver,digits=3)), Nse = $(round(Nse_Aver,digits=3))  ~")

			println("    	~  Rmse_θΨ = $(round(Rmse_θΨ_Aver,digits=4)),  Rlmse_KΨ = $(round(Rmse_KΨ_Aver,digits=4)), Rmse = $(round(Rmse_Aver,digits=4))  ~ \n")

	return hydro, hydroOther
	end  # function: HYPIXOPT_START


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OF_HYPIX
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function OF_HYDROLAB(X, hydro, iZ, K_KΨ, N_KΨ, N_θΨ, optim, optionHydro, θ_θΨ, Ψ_KΨ, Ψ_θΨ)

			# New optimized param which are put into the matching veg or hydro parameters
				hydro = hydrolabOpt.PARAM_2_hydro(hydro, iZ, optim, optionHydro, X)
		
			# Weighted Objective Function
				Of, Of_θΨ, Of_Kunsat = OfHydrolab.OF_WRC_KUNSAT(iZ, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro, optim, optionHydro) 
				
		return Of
		end  # function: OF_HYPIX


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PARAM
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function PARAM_2_hydro(hydro, iZ, optim, optionHydro, X)
		
		for iParam = 1:optim.NparamOpt
			# Determening if parameters are Log transformed
				if (optim.ParamOpt_LogTransform[iParam]) && !(optim.ParamOpt[iParam]=="Ψm" && optionHydro.σ_2_Ψm == :Constrained)
					Paramₐ = expm1(X[iParam])
				else
					Paramₐ = X[iParam]
				end  # if: optim.ParamOpt_LogTransform

			# Getting the current values of every layer of the hydro parameter of interest
				vectParam = getfield(hydro, Symbol(optim.ParamOpt[iParam]))

			# Updating the value of the parameters for the soil wanting to optimize by keeping the values constant
				vectParam[iZ] = Paramₐ

			# Putting the updated hydro param into hydro
				setfield!(hydro, Symbol(optim.ParamOpt[iParam]), vectParam)
		end # for loop


		# ==================== SPECIAL CASE ====================

		# RELATIONSHIP BETWEEN σ AND Ψm
		if (optionHydro.σ_2_Ψm ≠ :No) && ("Ψm" ∈ optim.ParamOpt)
			hydro = hydroRelation.FUNCTION_σ_2_Ψm_SOFTWARE(hydro, iZ, option.hydro; Pσ=3.0)
		end # option.hydro.σ_2_Ψm ≠ :No

		#  <>=<>=<>=<>=<>=<> Relationship between σ and θr
		if option.hydro.θrOpt==:σ_2_θr && ("θr" ∉ optim.ParamOpt) && ("σ" ∈ optim.ParamOpt)
			hydro.θr[iZ] = hydroRelation.σ_2_θr(hydro, iZ)
		end

		# Converting θsMacMat_ƞ -> θsMacMat
		if  optionHydro.HydroModel == :Kosugi
			hydro.θsMacMat[iZ] = hydro.θsMacMat_ƞ[iZ] * hydro.θs[iZ]
		end

	return hydro
	end  # function: PARAM

end  # module hypixOpt
# ............................................................