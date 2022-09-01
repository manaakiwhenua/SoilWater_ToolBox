
# =============================================================
#		MODULE: pond
# =============================================================
module pond
	import ..kunsat, ..param, ..cst
	export PONDING_SORPTIVITY!

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PONDING
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# function PONDING(discret, ΔHpond, hydro, iT::Int,  ΔPr, ΔT, Ψ,  θ)

		# 	ΔKs = ΔT[iT] * K_Aver[1]

		# 	ΔHpond[iT] = (ΔPr[iT] + ΔHpond[iT-1] - ΔKs * ((Ψ[iT,1] / discret.ΔZ_Aver[1]) + param.hyPix.Cosα) ) / ((ΔKs / discret.ΔZ_Aver[1]) + 1.0)  

		# 	ΔHpond[iT] = max(ΔHpond[iT], 0.0)

		# 	# More ponding will occure if the amont of infiltration is greater than the current storage
		# 	ΔHsurplus = max(ΔPr[iT] + ΔHpond[iT-1] - ΔHpond[iT] - discret.ΔZ[1] * max(hydro.θs[1] - θ[iT-1,1], 0.0) , 0.0)

		# 	ΔHpond[iT] += ΔHsurplus
			
		# 	return ΔHpond
		# end  # function: PONDING

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PONDING_SORPTIVITY!
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function PONDING_SORPTIVITY!(discret, hydro, iT, Sorptivity, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ)

				Bparam = (2.0 - cst.β) / 3.0 + (kunsat.Ψ_2_KUNSAT(Ψ[iT-1,1], 1, hydro) / hydro.Ks[1]) * (1.0 + cst.β) / 3.0
				
				Infilt_Max =  (Sorptivity * √ΔT[iT] + Bparam * hydro.Ks[1] * ΔT[iT]) * param.hyPix.Cosα

			# Reduction of infiltration rate to avoid that too much water infiltrates into layer 1
				Infilt_Max = min(discret.ΔZ[1] * (hydro.θs[1] - θ[iT-1,1]) + ΔSink[iT,1], Infilt_Max)

			# ΔHpond is computed as
				ΔHpond[iT] = max(ΔPr[iT] + ΔHpond[iT-1] - Infilt_Max, 0.0)
				
			return ΔHpond
		end  # function: PONDING
	
end  # module pond
# ............................................................