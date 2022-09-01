# =============================================================
#		MODULE: psdInitiate
# =============================================================
module psdInitialize
	import ..psdStruct, ..psdThetar, ..param
	export PSD_INITIALIZE
	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PSD_INITIALIZE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function PSD_INITIALIZE(∑Psd, hydro, hydroPsd, N_Psd, N_SoilSelect)			
		# Compute new N_Psd to take into account when ∑Psd 
		# Correction for N_Psd such that to determine the real maximum Rpart size
			for iZ=1:N_SoilSelect
				N_Psd_New = 1
				for iPsd = 1:N_Psd[iZ]
					if ∑Psd[iZ, iPsd] >= 0.99999
						N_Psd_New = iPsd
						break
					end
				end #  iPsd = 1:N_Psd[iZ]

				N_Psd[iZ] = N_Psd_New
			end # looping over soils

		# MAximum number of Psd data
			N_Psd_Max = maximum(N_Psd[1:N_SoilSelect])


		# Compute PSD from ∑PSD
			Psd = zeros(Float64, (N_SoilSelect, N_Psd_Max))
			for iZ=1:N_SoilSelect
				Psd[iZ, 1:N_Psd[iZ]] = ∑PSD_2_PSD(∑Psd[iZ,1:N_Psd[iZ]], N_Psd[iZ])
			end # iZ=1:N_SoilSelect

		# DERIVING THE STRUCTURE PARAMETERS
			paramPsd = psdStruct.PSDSTRUCT(N_SoilSelect)

		# COMPUTING θr FROM PSD DATA
		hydroPsd, paramPsd = psdThetar.PSD_2_θr(∑Psd, hydro, hydroPsd, N_SoilSelect, paramPsd)
			
		# COMPUTING θs FROM TOTAL POROSITY
			θs_Psd = Array{Float64}(undef, N_SoilSelect)
			for iZ=1:N_SoilSelect
				θs_Psd[iZ] = param.hydro.Coeff_Φ_2_θs * hydroPsd.Φ[iZ]
				hydroPsd.θs[iZ] = θs_Psd[iZ]
				hydroPsd.θr[iZ] = paramPsd.θr_Psd[iZ]
			end 
		
	return N_Psd, N_Psd_Max, Psd, θs_Psd, hydroPsd, paramPsd
	end  # function: PSD_INITIALIZE


	# =========================================
	#          ∑PSD -> PSD
	# =========================================
		function ∑PSD_2_PSD(∑Psd, N_Psd)
			Psd = zeros(Float64, N_Psd)

			Psd[1] = ∑Psd[1]
			for iRpart =2:N_Psd
				Psd[iRpart] = ∑Psd[iRpart] - ∑Psd[iRpart-1]
			end
		return Psd
		end # function ∑PSD_2_PSD
	
end  # module psdInitiate
# ............................................................