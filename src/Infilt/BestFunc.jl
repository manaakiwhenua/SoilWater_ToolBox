# =============================================================
#		MODULE: bestFunc
# =============================================================
module bestFunc
	import ..option, ..sorptivity, ..wrc, ..kunsat, ..option
	export  BEST_UNIVERSAL_START, CONVERT_3D_2_1D, CONVERT_1D_2_3D

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : BEST
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function BEST_UNIVERSAL_START(∑Infilt_3D, hydroInfilt, infiltOutput, infiltParam, iZ, N_Infilt, T; θ_Ini=infiltParam.θ_Ini[iZ])

			# Initializing
				Se_Ini = wrc.θ_2_Se(θ_Ini, iZ, hydroInfilt)

				Kr_θini = (kunsat.Se_2_KUNSAT(Se_Ini, iZ, hydroInfilt)) / hydroInfilt.Ks[iZ]

				Sorptivity = sorptivity.SORPTIVITY(θ_Ini, iZ, hydroInfilt)

				A = bestFunc.A(θ_Ini, hydroInfilt.θs[iZ], iZ, infiltParam)

				B = bestFunc.B(iZ, Kr_θini, infiltParam)

				T_TransSteady = bestFunc.TIME_TRANS_STEADY(B, hydroInfilt.Ks[iZ], Sorptivity)

				for iT = 1:N_Infilt[iZ]
					∑Infilt_3D[iZ, iT] = BEST_UNIVERSAL(iZ, A, B, Sorptivity, T[iZ,iT], T_TransSteady, hydroInfilt, infiltParam)
				end  # for iT=1:N_Infilt[iZ]

				return ∑Infilt_3D, T_TransSteady
		end # function: BEST_UNIVERSAL_START

		# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
		# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : BEST_UNIVERSAL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function BEST_UNIVERSAL(iZ, A, B, Sorptivity, T, T_TransSteady, hydroInfilt, infiltParam)

			if option.infilt.DataSingleDoubleRing == :Single #<>=<>=<>=<>=<>

				if T ≤ T_TransSteady
					return ∑Infilt_3D = bestFunc.INFILTRATION_3D_TRANSIT(A, B, hydroInfilt.Ks[iZ], Sorptivity, T)
				else
					return ∑Infilt_3D = bestFunc.INFILTRATION_3D_STEADY(A, B, iZ, hydroInfilt.Ks[iZ], Sorptivity, T, infiltParam, T_TransSteady)
				end # T <= T_TransSteady

			elseif option.infilt.DataSingleDoubleRing == :Double  #<>=<>=<>=<>=<>
				if T ≤ T_TransSteady
					return ∑Infilt_3D = bestFunc.INFILTRATION_1D_TRANSIT(B, hydroInfilt.Ks[iZ], Sorptivity, T)
				else
					return ∑Infilt_3D = bestFunc.INFILTRATION_1D_STEADY(B, iZ, hydroInfilt.Ks[iZ], Sorptivity, T, infiltParam, T_TransSteady)
				end # T <= T_TransSteady

			end # option.∑Infilt_3D.Dimension
			
		end  # function: BEST_UNIVERSAL


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION_3D_TRANSIT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATION_3D_TRANSIT(A, B, Ks, Sorptivity, T)
			return Sorptivity * (T ^ 0.5) + (A * (Sorptivity ^ 2.0) + B * Ks) * T
		end  # function: INFILTRATION_3D_TRANSIT
	

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION_1D_TRANSIT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATION_1D_TRANSIT(B, Ks, Sorptivity, T)
			return Sorptivity * √T + B * Ks * T
		end # function: INFILTRATION_1D_TRANSIT


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION_3D_STEADY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATION_3D_STEADY(A, B, iZ, Ks, Sorptivity, T, infiltParam, T_TransSteady)
			if option.infilt.bestFunc.Continous == true
				return bestFunc.INFILTRATION_3D_TRANSIT(A, B, Ks, Sorptivity, T_TransSteady)  + (A * (Sorptivity ^ 2.0) + Ks) * (T - T_TransSteady)

			elseif option.infilt.bestFunc.Continous  == false
				return (A * (Sorptivity ^ 2.0) + Ks) * T + bestFunc.C(B, infiltParam, iZ) * (Sorptivity ^ 2.0) / Ks
			end
		end  # function: INFILTRATION_3D_STEADY
# 

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION_1D_STEADY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTRATION_1D_STEADY(B, iZ, Ks, Sorptivity, T, infiltParam, T_TransSteady)
			if option.infilt.bestFunc.Continous == true
				return bestFunc.INFILTRATION_1D_TRANSIT(B, Ks, Sorptivity, T_TransSteady)  + Ks * (T - T_TransSteady)
				
			elseif option.infilt.bestFunc.Continous == false
				return Ks * T + bestFunc.C(B,  infiltParam, iZ) * (Sorptivity ^ 2.0) / Ks
			end
		end  # function: INFILTRATION_1D_STEADY


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : CONVERT_3D_2_1D
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function CONVERT_3D_2_1D(∑Infilt_3D, ∑Infilt_1D, hydroInfilt, infiltParam, iZ, N_Infilt, T; θ_Ini= infiltParam.θ_Ini[iZ])
				
			Sorptivity = sorptivity.SORPTIVITY(θ_Ini, iZ, hydroInfilt)

			A = bestFunc.A(θ_Ini, hydroInfilt.θs[iZ], iZ, infiltParam)

			for iT=1:N_Infilt[iZ]
				∑Infilt_1D[iZ,iT] = ∑Infilt_3D[iZ,iT] - A  * T[iZ,iT] * Sorptivity ^ 2.0
			end # iT

			return ∑Infilt_1D
		end  # function: CONVERT_3D_2_1D
		
		
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : CONVERT_1D_2_3D
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function CONVERT_1D_2_3D(∑Infilt_3D, ∑Infilt_1D, hydroInfilt, infiltParam, iZ, N_Infilt, T; θ_Ini= infiltParam.θ_Ini[iZ])

			Sorptivity = sorptivity.SORPTIVITY(θ_Ini, iZ, hydroInfilt)

			A = bestFunc.A(θ_Ini, hydroInfilt.θs[iZ], iZ, infiltParam)

			for iT=1:N_Infilt[iZ]
				∑Infilt_3D[iZ,iT] = ∑Infilt_1D[iZ,iT] + A * (Sorptivity ^ 2.0) * T[iZ,iT]
			end # iT

			return ∑Infilt_3D
		end  # function: CONVERT_1D_2_3D


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : TIME_TRANS_STEADY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function TIME_TRANS_STEADY(B, Ks, Sorptivity)
			return ( Sorptivity / (Ks * 2.0 * (1.0 - B)) ) ^ 2.0
		end # function: TIME_TRANS_STEADY


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : A
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function A(θ_Ini, θs, iZ, infiltParam)
			return  infiltParam.γ[iZ] / ( infiltParam.RingRadius[iZ] * (θs - θ_Ini)) # Units [mm-1]
		end  # function: A
	

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : B
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function B(iZ, Kr_θini,  infiltParam)
			return (2.0 -  infiltParam.β[iZ]) / 3.0 + Kr_θini * (1.0 + infiltParam.β[iZ]) / 3.0
		end # function: B


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : C
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function C(B, infiltParam, iZ)
			return log(1.0 /  infiltParam.β[iZ]) * (1.0 +  infiltParam.β[iZ]) / (6.0 * (1.0 -  infiltParam.β[iZ]) * (1.0 - B) )
		end # function: C

end # MODULE: bestFunc