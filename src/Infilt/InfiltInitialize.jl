# =============================================================
#		MODULE: HYDRO INITIALIZE
# =============================================================
module infiltInitialize

	import ..hydroStruct, ..psdThetar, ..param, ..timeTransSteady, ..infiltStruct, ..option
	export INFILT_INITIALIZE

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILT_INITIALIZE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILT_INITIALIZE(∑Infilt_Obs, ∑Psd, hydroInfilt, infiltParam, N_Infilt, N_SoilSelect, Tinfilt)

			# CORRECTION FOR θr & θs
				for iZ=1:N_SoilSelect
					# θr computation
						# if option.globalopt.Psd
							hydroInfilt.θr[iZ] = psdThetar.PSD_2_θr_FUNC(∑Psd, hydroInfilt, iZ)
						# else
						# 	hydroInfilt.θr[iZ] = param.hydro.θr
						# end # option.globalopt.Psd

					# θr < θ_Ini
						infiltParam.θ_Ini[iZ] = max(hydroInfilt.θr[iZ] + eps(), infiltParam.θ_Ini[iZ])

					# θs computation
						hydroInfilt.θs[iZ] = param.hydro.Coeff_Φ_2_θs * hydroInfilt.Φ[iZ]

					# Checking for θ_Ini
						infiltParam.θ_Ini[iZ] = min(infiltParam.θ_Ini[iZ], hydroInfilt.θs[iZ] - 0.0015) 
				end  # for iZ=1:N_SoilSelect

			# TIME FLUX CORRECTION
				N_Infilt_Max = maximum(N_Infilt[1:N_SoilSelect])

				T = Array{Float64}(undef, (N_SoilSelect, N_Infilt_Max))
				for iZ=1:N_SoilSelect
					# T[iZ,1] = ( (Tinfilt[iZ,1] ^0.5 + (0.0)^0.5) / 2.0 ) ^ 2.0
					T[iZ,1] = Tinfilt[iZ,1]

					for iInfilt=2:N_Infilt[iZ]
						T[iZ,iInfilt] = ( (Tinfilt[iZ,iInfilt] ^ 0.5 + Tinfilt[iZ,iInfilt-1] ^ 0.5) / 2.0 ) ^ 2.0
					end	
				end #  iZ=1:N_SoilSelect

			# STRUCTURE OF INFILTRATION PARAMETERS
				infiltOutput = infiltStruct.INFILTSTRUCT(N_SoilSelect)

			# DETERMENING WHEN STEADY STATE OCCURES
				infiltOutput = timeTransSteady.∑INFIlT_2_TIMETRANSSTEADY(T, N_SoilSelect, N_Infilt, infiltOutput, ∑Infilt_Obs) 

			# Initializing Infilt		
				∑Infilt_3D = Array{Float64}(undef, (N_SoilSelect, N_Infilt_Max))
				∑Infilt_1D = Array{Float64}(undef, (N_SoilSelect, N_Infilt_Max))

			return T, infiltOutput, hydroInfilt, ∑Infilt_3D, ∑Infilt_1D

		end  # function: INFILT_INITIALIZE


end # module hydroInitialize
# ............................................................