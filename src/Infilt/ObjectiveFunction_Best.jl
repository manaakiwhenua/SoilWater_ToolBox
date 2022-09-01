# =============================================================
#		MODULE: bestOpt
# =============================================================
module ofBest
	import ..stats, ..bestFunc
	export OF_BEST

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OF_BEST
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function OF_BEST(∑Infilt_3D, ∑Infilt_Obs, hydroInfilt, infiltOutput, infiltParam, iZ, N_Infilt, T; W=0.3)

		∑Infilt_3D, T_TransSteady = bestFunc.BEST_UNIVERSAL_START(∑Infilt_3D, hydroInfilt, infiltOutput, infiltParam, iZ, N_Infilt, T)

		iT_TransSteady = infiltOutput.iT_TransSteady_Data[iZ]

		Nse_Trans = stats.NASH_SUTCLIFE_MINIMIZE( ∑Infilt_Obs[iZ,1:iT_TransSteady], ∑Infilt_3D[iZ, 1:iT_TransSteady]; Power=2.0)

		Nse_Steady = stats.NASH_SUTCLIFE_MINIMIZE( log10.(∑Infilt_Obs[iZ,iT_TransSteady+1:N_Infilt[iZ]]), log10.(∑Infilt_3D[iZ,iT_TransSteady+1:N_Infilt[iZ]]); Power=2.0)

		Penalty = abs(∑Infilt_Obs[iZ,N_Infilt[iZ]] - ∑Infilt_3D[iZ,N_Infilt[iZ] ]) /  ∑Infilt_Obs[iZ,N_Infilt[iZ]]

		return Nse = W * Nse_Trans + (1.0 - W) * Nse_Steady  + Penalty
		
	end  # function: OF_BEST


	
end  # module: bestOpt
# ............................................................