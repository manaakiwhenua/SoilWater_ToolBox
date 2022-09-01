# =============================================================
#		MODULE: hydroRealation
# =============================================================
module hydroRelation
using BlackBoxOptim
import ..param, ..tool
export σ_2_Ψm, σ_2_θr

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : σ_2_Ψm(iZ, hydro)
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function σ_2_Ψm(σ₁, Ψσ, Ψm_Min, Ψm_Max; Pσ=3.0)
			Ψm = Ψσ * exp(σ₁*Pσ)
			return Ψm = max(min(Ψm , Ψm_Max), Ψm_Min)
		end # function: σ_2_Ψm


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : σ_2_θr
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function FUNCTION_σ_2_Ψm_SOFTWARE(hydro₂, iZ, option₂; Pσ=3.0)
			if (option₂.σ_2_Ψm == :Constrained)
				Ψm_Min = hydroRelation.σ_2_Ψm(hydro₂.σ[iZ], param.hydro.kg.Ψσ_Min, hydro₂.Ψm_Min[iZ], hydro₂.Ψm_Max[iZ])

				Ψm_Max = hydroRelation.σ_2_Ψm(hydro₂.σ[iZ], param.hydro.kg.Ψσ_Max, hydro₂.Ψm_Min[iZ], hydro₂.Ψm_Max[iZ])
				
				hydro₂.Ψm[iZ] = tool.norm.∇NORM_2_PARAMETER(hydro₂.Ψm[iZ], Ψm_Min, Ψm_Max)

			elseif (option₂.σ_2_Ψm == :UniqueRelationship) # <>=<>=<>=<>=<>
				hydro₂.Ψm[iZ] = hydroRelation.σ_2_Ψm(hydro₂.σ[iZ], param.hydro.Ψσ, hydro₂.Ψm_Min[iZ], hydro₂.Ψm_Max[iZ])

			end #option.infilt.σ_2_Ψm

		return hydro₂
		end


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : σ_2_θr
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function σ_2_θr(hydro₂, iZ; α₁=17.5, α₂=4.0)
			σ_η = (hydro₂.σ[iZ] - hydro₂.σ_Min[iZ]) / (hydro₂.σ_Max[iZ] - hydro₂.σ_Min[iZ]) 
			
			return θr =  (hydro₂.θr_Max[iZ] * (1.0 - exp(-α₁ * σ_η ^ α₂))) / (1.0 - exp(-α₁ * 1.0 ^ α₂))
		end  # function: σ_2_θr


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : σ_2_Ψm_Old(iZ, hydro)
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# function σ_2_Ψm(iZ::Int, hydro; σ_Min= param.hydro.kg.σ_Min, σ_Max=param.hydro.kg.σ_Max )
		# 	hydro.σ[iZ] = param.hydro.kg.Pσ_1 * ( log(hydro.Ψm[iZ]) - 1.0 ) ^ param.hydro.kg.Pσ_2
		# 	hydro.σ[iZ] = max(min(hydro.σ[iZ], σ_Max), σ_Min)
		# 	return hydro
		# end  # function: σ_2_Ψm


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : σ_2_Ψm_OLD(iZ, hydro)
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# function σ_2_Ψm(iZ::Int, hydro)
			# hydro.Ψm[iZ] = exp(exp( inv(param.hydro.kg.Pσ_2) * log(hydro.σ[iZ] / param.hydro.kg.Pσ_1)) + 1.0) #[mm]
		# 	hydro.Ψm[iZ] = max(min(hydro.Ψm[iZ] , param.hydro.kg.Ψm_Max), param.hydro.kg.Ψm_Min)
		# 	return hydro
		# end # function: σ_2_Ψm


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OPTIMIZATION_σ_2_Ψm_oLD
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# function OPTIMIZATION_σ_2_Ψm(Ψm_Obs, σ_Obs)

		# 	function OF_σ_2_Ψm(Ψm_Obs, σ_Obs, Pσ_1, Pσ_2)
		# 		Of = 0.0
		# 		i = 1
		# 		for σ in σ_Obs
		# 			Ψm_Sim = exp( exp(inv(Pσ_2) * log(σ / Pσ_1)) + 1.0 ) #[mm]
		# 			Ψm_Sim = max(min(Ψm_Sim , param.hydro.kg.Ψm_Max), param.hydro.kg.Ψm_Min)
		# 			Of += (log(Ψm_Sim) - log(Ψm_Obs[i]))^2.  
		# 			i += 1
		# 		end # for σ in σ_Obs
		# 		return Of 
		# 	end # function: OF_σ_2_Ψm

		# 	Optimization = BlackBoxOptim.bboptimize(Pσ ->  OF_σ_2_Ψm(Ψm_Obs, σ_Obs, Pσ[1], Pσ[2]); SearchRange = [(0., 100.), (0., 100.)], NumDimensions=2, TraceMode=:silent)
				
		# 	Pσ_1 = BlackBoxOptim.best_candidate(Optimization)[1]
		# 	Pσ_2 = BlackBoxOptim.best_candidate(Optimization)[2]

		# 	println("Pσ_1 = $Pσ_1, Pσ_2 = $Pσ_2")
		# 	return Pσ_1, Pσ_2
		# end # function OPTIMIZATION_σ_2_Ψm

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : σ_2_Ψm(iZ, hydro)
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# function σ_2_Ψm(σ₁; Ψσ=param.hydro.kg.Ψσ, Pσ=1.0, Ψm_Min=param.hydro.kg.Ψm_Min, Ψm_Max=param.hydro.kg.Ψm_Max)
		# 	Ψm = Ψσ * exp(σ₁ * (σ₁ + Pσ))
		# 	 return Ψm = max(min(Ψm , Ψm_Max), Ψm_Min)
		# end # function: σ_2_Ψm
	
end  # module: hydroRealation
# ............................................................