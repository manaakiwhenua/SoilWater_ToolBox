# =============================================================
#		module: ΨminΨmax
# =============================================================
module ΨminΨmax

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ΨMINΨMAX
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function ΨMINΨMAX(θs, θsMacMat, σ, σMac, Ψm, ΨmMac; Pσ=4.0)

		# Ψ_Min
			# Have we got macropore ?
			if θs - θsMacMat > 0.001
				Ψ_Min = exp( log(ΨmMac) - σMac * Pσ) 
			else
				Ψ_Min = exp( log(Ψm) - σ * Pσ)
			end  # if: θs -θsMacMat > 0.01

		#  Ψ_Max
			Ψ_Max = exp(log(Ψm) + σ * Pσ)

		return Ψ_Max, Ψ_Min
	end  # function: ΨMINΨMAX
	
end  # module ΨminΨmax
# ............................................................