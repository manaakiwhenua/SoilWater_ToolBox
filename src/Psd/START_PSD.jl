module psd
	import ..option, ..param, ..stats, ..psdInitialize, ..psdOpt, ..psdFunc

	# ======================================================================================
	#          PSD_START Slow initialization
	# ======================================================================================
	function START_PSD(∑Psd, hydro, hydroPsd, N_Psd, N_SoilSelect, N_θΨ, Rpart, θ_θΨ, Ψ_θΨ)

		if  option.psd.OptimizePsd == :OptAllSoil && option.globalopt.θΨ == :No
			error("PSD error: if  option.psd.OptimizePsd == :OptAllSoil && option.globalopt.θΨ ≠ :No ")
		end 

		# INITIATING THE PSD DATA		
		N_Psd, N_Psd_Max, Psd, θs_Psd, hydroPsd, paramPsd = psdInitialize.PSD_INITIALIZE(∑Psd, hydro, hydroPsd, N_Psd, N_SoilSelect)

		if option.psd.OptimizePsd == :Run  # <>=<>=<>=<>=<>=<>=<>=<>=<>=<> 
			θ_Rpart, Ψ_Rpart = psdOpt.PSD_RUN_ALLMODEL(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, paramPsd.θr_Psd, paramPsd, hydroPsd)
			
			if option.globalopt.θΨ ≠ :No # <>=<>=<>=<>=<>=<>=<>=<>=<>=<> 
				paramPsd.Nse, Nse_Mean_Run, Nse_Std_Run = stats.NASH_SUTCLIFFE_θΨ(N_SoilSelect, N_Psd, Ψ_Rpart, θ_Rpart, hydroPsd)
				println("\n    == RUN the PSD parameters with prescribed parameters ==")
				println("    ~ Nse_Mean_Run=$Nse_Mean_Run, Nse_Std_Run=$Nse_Std_Run ~")
			end
		
			# && option.globalopt.θΨ ≠ :No
		elseif option.psd.OptimizePsd == :OptSingleSoil # <>=<>=<>=<>=<>=<>
			if option.psd.Model == :IMP # <>=<>=<>=<>=<>=<>=<>=<>=<>=<> 
				paramPsd, θ_Rpart, Ψ_Rpart, Nse_Mean_SingleOpt, Nse_Std_SingleOpt = psdOpt.imp.OPTIMIZATION_SINGLE_SOIL(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, paramPsd.θr_Psd, paramPsd, hydro)

			elseif option.psd.Model == :Chang2019Model # <>=<>=<>=<>=<>=<>=<>=<>=<>=<> 
				paramPsd, θ_Rpart, Ψ_Rpart, Nse_Mean_SingleOpt, Nse_Std_SingleOpt = psdOpt.chang.OPTIMIZATION_SINGLE_SOIL(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, paramPsd.θr_Psd, paramPsd, hydro)
			end

			println("\n    == Optimizing the PSD parameters individually for each soil ==")
			println("    	~ Nse_Mean_SingleOpt=$Nse_Mean_SingleOpt,  Nse_Std_SingleOpt=$Nse_Std_SingleOpt ~ \n")

			# && option.globalopt.θΨ ≠ :No
		elseif option.psd.OptimizePsd == :OptAllSoil # <>=<>=<>=<>=<>=<>
			if option.psd.Model == :IMP # <>=<>=<>=<>=<>=<>=<>=<>=<>=<> 
				paramPsd, θ_Rpart, Ψ_Rpart, Nse_Mean_OptAllSoil, Nse_Std_OptAllSoil = psdOpt.imp.OPTIMIZATION_ALL_SOIL(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, paramPsd.θr_Psd, paramPsd, hydro)

			elseif option.psd.Model == :Chang2019Model # <>=<>=<>=<>=<>=<>=<>=<>=<>=<> 
				paramPsd, θ_Rpart, Ψ_Rpart, Nse_Mean_OptAllSoil, Nse_Std_OptAllSoil = psdOpt.chang.OPTIMIZATION_ALL_SOIL(N_Psd_Max, N_SoilSelect, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, paramPsd.θr_Psd, paramPsd, hydro)
			end

			println("\n    == Optimizing the PSD parameters for all soils ==")
			println("    	~ Nse_Mean_OptAllSoil=$Nse_Mean_OptAllSoil,  Nse_Std_OptAllSoil=$Nse_Std_OptAllSoil ~")

  		else
     		error("  $(option.psd.OptimizePsd) not found ")
		end # option.psd.OptimizePsd

	return paramPsd, N_Psd, θ_Rpart, Ψ_Rpart, Psd, hydroPsd	  	
	end # function START_PSD

end # module PSD