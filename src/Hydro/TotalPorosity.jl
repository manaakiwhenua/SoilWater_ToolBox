# =============================================================
#		MODULE: totalPorosity
# =============================================================
module Φ
	export  ρB_2_Φ

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Φ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ρB_2_Φ(N_SoilSelect, RockW, ρ_Rock, ρbSoil, ρp_Fine)

			Φ = Array{Float64}(undef, (N_SoilSelect))

			for iZ=1:N_SoilSelect
				# Vrock = RockW[iZ] / ρ_Rock[iZ]
				Φ[iZ] = 1.0 - (RockW[iZ] * ρbSoil[iZ] / ρ_Rock[iZ]) - ((1.0 - RockW[iZ]) *ρbSoil[iZ] / ρp_Fine[iZ])
			end # for
			
			return Φ
		end  # function: Φ
	
end  # module: totalPorosity
# ............................................................