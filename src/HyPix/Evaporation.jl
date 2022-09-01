# =============================================================
#		MODULE: evaporation
# =============================================================
module evapo
	import ..wrc
	export N_IEVAPO, EVAPORATION!

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :  COMPUTTING N OF EVAPORATION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function N_IEVAPO(N_iZ::Int64, veg, Z)::Int64
			iZ = 1
			while iZ ≤ N_iZ && veg.Zevapo ≥ Z[iZ]
				iZ += 1
			end
			return N_iEvapo = iZ - 1
		end  # function EVAPORATION


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : EVAPORATION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function EVAPORATION!(hydro, iT, ΔEvaporation, ΔPet_Evap, θ)
			# Se_Max = 0.0
			# iZ_Evapo = 1

			# for iZ = 1:N_iEvapo
			# 	Se = wrc.θ_2_Se(θ[iT-1,iZ], iZ, hydro)
			# 	if Se > Se_Max
			# 		Se_Max = Se
			# 		iZ_Evapo = iZ
			# 	end
			# end

			# iZ_Evapo = 1
			
			ΔEvaporation[iT] = ΔPet_Evap * wrc.θ_2_Se(θ[iT-1,1], 1, hydro)
			
			return ΔEvaporation
		end  # function EVAPORATION
	
end  # module evapo
# ............................................................