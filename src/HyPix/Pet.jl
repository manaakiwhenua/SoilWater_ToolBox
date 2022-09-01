# =============================================================
#		MODULE: et
# =============================================================
module pet
	import ..option
	export BEER_LAMBERT_LAW

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :  BEER_LAMBERT_LAW
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function BEER_LAMBERT_LAW(iT, Lai, ΔPet, veg)

		ΔPet_Evap = ΔPet[iT] * exp(-Lai * veg.ExtinctCoefRadiation)
		 
		ΔPet_Transp = ΔPet[iT] - ΔPet_Evap

		return ΔPet_Evap, ΔPet_Transp
	end  # function BEER_LAMBERT_LAW
	
end  # module et
# ............................................................