module kunsat
	import ..option, ..wrc
	export Ψ_2_KUNSAT, Se_2_KUNSAT, θ_2_KUNSAT

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Ψ_2_KUNSAT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function Ψ_2_KUNSAT(Ψ₁, iZ::Int64, hydroParam)
			if option.hydro.HydroModel == :Kosugi
				return Kunsat = kunsat.kg.Ψ_2_KUNSAT(Ψ₁, iZ::Int64, hydroParam)
			elseif option.hydro.HydroModel == :Vangenuchten || option.hydro.HydroModel == :VangenuchtenJules
				return Kunsat = kunsat.vg.Ψ_2_KUNSAT(Ψ₁, iZ::Int64, hydroParam)
			elseif option.hydro.HydroModel == :BrooksCorey
				return Kunsat = kunsat.bc.Ψ_2_KUNSAT(Ψ₁, iZ::Int64, hydroParam)
			elseif option.hydro.HydroModel == :ClappHornberger
				return Kunsat = kunsat.ch.Ψ_2_KUNSAT(Ψ₁, iZ::Int64, hydroParam)
			else
				error("$(option.hydro.HydroModel) model for Ψ_2_KUNSAT is not yet available")
			end
		end # function Ψ_2_KUNSAT


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ΨSE_2_KUNSAT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  function ΨSE_2_KUNSAT(Ψ₁, Se, iZ::Int64, hydroParam)
			if option.hydro.HydroModel == :Kosugi
				return Kunsat = kunsat.kg.ΨSE_2_KUNSAT(Ψ₁, Se, iZ::Int64, hydroParam)
			else
				error("$(option.hydro.HydroModel) model for ΨSE_2_KUNSAT is not yet available")
			end
		end # function Ψ_2_KUNSAT


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Se_2_KUNSAT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  function Se_2_KUNSAT(Se, iZ::Int64, hydroParam)
			Se = max(min(Se, 1.0), 0.0)

			if option.hydro.HydroModel == :Kosugi
				return Kunsat = kunsat.kg.Se_2_KUNSAT(Se, iZ::Int64, hydroParam)
			elseif option.hydro.HydroModel == :Vangenuchten
				return Kunsat = kunsat.vg.Se_2_KUNSAT(Se, iZ::Int64, hydroParam)
			elseif option.hydro.HydroModel == :BrooksCorey
				return Kunsat = kunsat.bc.Se_2_KUNSAT(Se, iZ::Int64, hydroParam)
			elseif option.hydro.HydroModel == :ClappHornberger
				return Kunsat = kunsat.ch.Se_2_KUNSAT(Se, iZ::Int64, hydroParam)
			else
				error("$(option.hydro.HydroModel) model for Se_2_KUNSAT is not yet available")
			end
		end # function Se_2_KUNSAT
	

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Se_2_KUNSAT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  function Se_2_KR(Se, iZ::Int64, hydroParam)
			Se = max(min(Se, 1.0), 0.0)

			if option.hydro.HydroModel == :Kosugi
				return Kunsat = kunsat.kg.Se_2_KR(Se, iZ::Int64, hydroParam)
			else
				error("$(option.hydro.HydroModel) model for Se_2_KR is not yet available")
			end
		end # function Se_2_KUNSAT


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : θ_2_KUNSAT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  function θ_2_KUNSAT(θ₁, iZ::Int64, hydroParam)
			if option.hydro.HydroModel == :Kosugi
				Se = wrc.θ_2_Se(θ₁, iZ, hydroParam)
				return Kunsat = Se_2_KUNSAT(Se, iZ, hydroParam)
			else
				error("$(option.hydro.HydroModel) model for θ_2_KUNSAT is not yet available")
			end
		end # function θ_2_KUNSAT


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ∂K∂Ψ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  function ∂K∂Ψ(Ψ₁, iZ::Int64, hydroParam)
			if option.hydro.HydroModel == :Kosugi
				return ∂Kunsat = kunsat.kg.∂K∂Ψ(Ψ₁, iZ::Int64, hydroParam)
			elseif option.hydro.HydroModel == :Vangenuchten
				return ∂Kunsat = kunsat.vg.∂K∂Ψ(Ψ₁, iZ::Int64, hydroParam)
			elseif option.hydro.HydroModel == :BrooksCorey
				return ∂Kunsat = kunsat.bc.∂K∂Ψ(Ψ₁, iZ::Int64, hydroParam)
			elseif option.hydro.HydroModel == :ClappHornberger
				return ∂Kunsat = kunsat.ch.∂K∂Ψ(Ψ₁, iZ::Int64, hydroParam)
			else
				error("$(option.hydro.HydroModel) model for ∂K∂Ψ is not yet available")
			end
		end # function ∂K∂Ψ


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ∂K∂Ψ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  function ∂K∂θ(Ψ₁, iZ::Int64, hydroParam)
			if option.hydro.HydroModel == :Kosugi
				return ∂Kunsat = kunsat.kg.∂K∂θ(Ψ₁, iZ::Int64, hydroParam)
			else
				error("$(option.hydro.HydroModel) model for ∂K∂θ is not yet available")
			end
		end

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : θψ_2_K
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function θΨ_2_KUNSAT(Se_Max, iZ::Int64, hydroParam, RockFragment::Float64; TopsoilSubsoil="Topsoil")
				if option.hydro.HydroModel == :Kosugi
					return Kunsat = kunsat.kg.θΨ_2_KUNSAT(Se_Max, iZ::Int64, hydroParam, RockFragment::Float64; TopsoilSubsoil="Topsoil")
				else
					error("$(option.hydro.HydroModel) model for θψ_2_KUNSAT is not yet available")
				end

			end  # function: θψ_2_K

	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>

	# =============================================================
	#		MODULE KOSUGI
	# =============================================================
	module kg
		import ..option, ..wrc, ...cst, ...param
		import ForwardDiff, QuadGK
		import SpecialFunctions: erfc, erfcinv
		export Ψ_2_KUNSAT, Se_2_KUNSAT, ∂K∂Ψ, θΨ_2_KUNSAT

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : Ψ_2_KUNSAT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function Ψ_2_KUNSAT(Ψ₁, iZ::Int64, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ])

				Se = wrc.Ψ_2_SeDual(Ψ₁, iZ, hydroParam)

				KsMat = Ks * (θsMacMat - θr) / (θs - θr)			
				Kunsat_Mat =  KsMat * √Se * (0.5 * erfc(((log(Ψ₁ / Ψm)) / σ + σ) / √2.0)) ^ 2

				KsMac = Ks * (θs - θsMacMat) / (θs - θr)
				Kunsat_Mac =  KsMac * √Se * (0.5 * erfc(((log(Ψ₁ / ΨmMac)) / σMac + σMac) / √2.0)) ^ 2
				
				return Kunsat = Kunsat_Mat + Kunsat_Mac
			end # function Ψ_2_KUNSAT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : Ψ_2_KUNSAT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ΨSE_2_KUNSAT(Ψ₁, Se, iZ::Int64, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ])

				KsMat = Ks * (θsMacMat - θr) / (θs - θr)			
				Kunsat_Mat =  KsMat * √Se * (0.5 * erfc(((log(Ψ₁ / Ψm)) / σ + σ) / √2.0)) ^ 2

				KsMac = Ks * (θs - θsMacMat) / (θs - θr)
				Kunsat_Mac =  KsMac * √Se * (0.5 * erfc(((log(Ψ₁ / ΨmMac)) / σMac + σMac) / √2.0)) ^ 2
				
				return Kunsat = Kunsat_Mat + Kunsat_Mac
			end # function ΨSE_2_KUNSAT

		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : Se_2_KUNSAT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function Se_2_KUNSAT(Se, iZ::Int64, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ])

				Se = max( min(Se, 1.0), 0.0)

				KsMat = Ks * (θsMacMat - θr) / (θs - θr)
				Kunsat_Mat = KsMat * √Se * (0.5 * erfc( erfcinv(2.0 * Se) + σ / √2.0 )) ^ 2

				KsMac = Ks * (θs - θsMacMat) / (θs - θr)
				Kunsat_Mac = KsMac * √Se * (0.5 * erfc( erfcinv(2.0 * Se) + σMac / √2.0 )) ^ 2

				return Kunsat = Kunsat_Mat + Kunsat_Mac
			end # function Se_2_KUNSAT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : Se_2_Kr
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function Se_2_KR(Se, iZ::Int64, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ])

				Se = max( min(Se, 1.0), 0.0)

				KrMat = (θsMacMat - θr) / (θs - θr)
				Kr_Mat = KrMat * √Se * (0.5*erfc( erfcinv(2.0*Se) + σ / √2.0 )) ^ 2.0

				KrMac = (θs - θsMacMat) / (θs - θr)
				Kr_Mac = KrMac * √Se * (0.5* erfc( erfcinv(2.0*Se) + σMac / √2.0 ))^2.0

				return Kr = Kr_Mat + Kr_Mac
			end # function: Se_2_KR


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂K∂Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂K∂Ψ(Ψ₁, iZ::Int64, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ])
				
				ψ =fill(0.0::Float64, 1) 

				∂K∂Ψ_Numerical(ψ::Vector) = Ψ_2_KUNSAT(abs(ψ[1]), iZ, hydroParam)
				
				ψ[1] = Ψ₁
				
				Func_∂K∂Ψ_Numerical = ψ -> ForwardDiff.gradient(∂K∂Ψ_Numerical, ψ)			
				∂K∂Ψ = Func_∂K∂Ψ_Numerical(ψ)[1]
				
				if isnan(∂K∂Ψ)
					∂K∂Ψ = 0.0
				end
				return ∂K∂Ψ 
			end # function: ∂K∂Ψ

			
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : ∂K∂θ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂K∂θ(Ψ₁, Se, iZ::Int64, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ]) 

				P1 = 1.0 / sqrt(2.0)
				P2 = 0.125

				KsMat = Ks * (θsMacMat - θr) / (θs - θr)

				∂Kunsat_Mat∂θ = KsMat * √Se * exp(-(σ / √2.0 + erfcinv( 2.0*Se )) ^ 2 + erfcinv( 2.0*Se )^2)*erfc(σ/ √2.0 + erfcinv(2.0*Se)) + P2*KsMat*erfc(σ / √2.0 + erfcinv(2.0*Se)) ^2 / √Se

				KsMac = Ks * (θs - θsMacMat) / (θs - θr)

				∂Kunsat_Mac∂θ = KsMac*sqrt(Se)*exp(-(σMac/ √2.0 + erfcinv(2.0*Se)) ^ 2 + erfcinv(2.0*Se) ^2)*erfc(σMac / √2.0 + erfcinv(2.0*Se)) + P2*KsMac*erfc(P1*σMac + erfcinv(2.0*Se)) ^ 2 / √Se

			return ∂Kunsat_Mat∂θ + ∂Kunsat_Mac∂θ
			end #  ∂K∂θ


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : θψ_2_K
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		"""Pollacco, J.A.P., Webb, T., McNeill, S., Hu, W., Carrick, S., Hewitt, A., Lilburne, L., 2017. Saturated hydraulic conductivity model computed from bimodal water retention curves for a range of New Zealand soils. Hydrol. Earth Syst. Sci. 21, 2725–2737. https://doi.org/10.5194/hess-21-2725-2017"""

			function θΨ_2_KUNSAT(Se_Max, iZ::Int64, hydroParam, RockFragment::Float64; TopsoilSubsoil="Topsoil", θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψm=hydroParam.Ψm[iZ], σ=hydroParam.σ[iZ], θsMacMat=hydroParam.θsMacMat[iZ], ΨmMac=hydroParam.ΨmMac[iZ], σMac=hydroParam.σMac[iZ], Ks=hydroParam.Ks[iZ], Rtol= 10^-8.0)

				# So it will be independent of the Rock Fragments
					if θs / (1.0 - RockFragment) - θsMacMat / (1.0 - RockFragment)  ≥ param.hydro.θs_θsMacMat # Bimodal
						FlagBimodal = true
					else
						FlagBimodal = false
					end

				# Correct for rock treshold	
				RockFragment_Treshold = 0.4
				if RockFragment > RockFragment_Treshold

					RockFragment2 = max(2 * RockFragment_Treshold - RockFragment, 0.0)

					θs = (θs / (1.0 - RockFragment)) * (1.0 - RockFragment2)
					
					θsMacMat = (θsMacMat / (1.0 - RockFragment)) * (1.0 - RockFragment2)

					θr = (θr / (1.0 - RockFragment)) * (1.0 - RockFragment2)
				end				

				if (TopsoilSubsoil=="Topsoil")
					if FlagBimodal
                  τ₁          = 5.007
                  τ₂          = 0.969
                  τ₃          = 0.787
                  τ₁Mac       = 4.734
                  τ₂Mac       = 0.511
						τ₃Mac       = 0.041
						σMac = 0.322
					else
                  τ₁ = 5.859
                  τ₂ = 0.967
                  τ₃ = 0.530
					end

				elseif (TopsoilSubsoil=="Subsoil")
					if FlagBimodal
                  τ₁          = 6.444
                  τ₂          = 0.859
                  τ₃          = 0.408
                  τ₁Mac       = 3.973
						τ₂Mac       = 0.642
                  τ₃Mac       = 0.729
						σMac = 1.272
					else # Unimodal  
                  τ₁          = 6.484
                  τ₂          = 0.854
                  τ₃          = 0.316		
					end
				end # TopsoilSubsoil

            # τ₁    = max(min(4.57 * (θsMacMat - θr) + 2.012, 6.4926), 1.4793)
            # τ₂    = max(min(0.050 * τ₁ + 0.242, 0.9), 0.1119)
            # τ₃    = max(min(-0.1479 * τ₁ + 1.011, 0.8), 0.0116)
            # τ₁Mac = max(min(-1.135 * (θsMacMat - θr) + 1.431, 4.9749), 0.0049)
            # τ₂Mac = max(min(0.405 * τ₂, 0.5896), 0.0122)
            # τ₃Mac = max(min(0.525 * τ₃, 0.7305), 0.0071)

				T1 = 10.0 ^ - τ₁
				T2 = 2.0 * (1.0 - τ₂)
				T3 =  1.0 / (1.0 - τ₃)

				Kunsat_Uni(Se) =  T1 * ((θsMacMat - θr) ^ T3) * ((cst.Y / Ψm) / (exp( erfcinv(2.0 * Se) * σ * √2.0 )) ) ^ T2
	
				if FlagBimodal
					T1Mac = 10.0 ^ - τ₁Mac 
					T2mac = 2.0 * (1.0 - τ₂Mac)
					T3mac = 1.0 / (1.0 - τ₃Mac)
					Kunsat_Bim(Se) = T1Mac * ((θs - θsMacMat) ^ T3mac) * ((cst.Y / ΨmMac) / ( exp( erfcinv(2.0 * Se) * σMac * √2.0))) ^ T2mac 
					
					Kunsat_Bimodal(Se) = Kunsat_Uni(Se) + Kunsat_Bim(Se) 
					return Kunsat_Model = cst.KunsatModel * QuadGK.quadgk(Se -> Kunsat_Bimodal(Se), 0.0, Se_Max)[1]
				else
					return Kunsat_Model = cst.KunsatModel * QuadGK.quadgk(Se -> Kunsat_Uni(Se), 0.0, Se_Max)[1]
				end   

				Kunsat_Model = min(max(hydroParam.Ks_Min[iZ], Kunsat_Model), hydroParam.Ks_Max[iZ])

			end # function θΨ_2_KUNSAT
	end # module kg 
	# =============================================================


	# =============================================================
	#		MODULE VAN GENUCHTEN
	# =============================================================
	module vg
		import ..option, ..wrc
		export Ψ_2_KUNSAT, Se_2_KUNSAT, ∂K∂Ψ

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION vg : Ψ_2_KUNSAT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function  Ψ_2_KUNSAT(Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψvg=hydroParam.Ψvg[iZ], N=hydroParam.N[iZ], Ks=hydroParam.Ks[iZ], Km=hydroParam.Km[iZ], L=0.5)
				M = 1.0 - Km / N

				Se = wrc.Ψ_2_SeDual(Ψ₁, iZ, hydroParam)
				return Kunsat = Ks * (Se^L) * ( 1.0 - (1.0 - Se ^ (1.0 / M) ) ^ M ) ^ 2.0
			end #function Ψ_2_KUNSAT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION vg : Se_2_KUNSAT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function  Se_2_KUNSAT(Se, iZ, hydroParam, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψvg=hydroParam.Ψvg[iZ], N=hydroParam.N[iZ], Ks=hydroParam.Ks[iZ], Km=hydroParam.Km[iZ], L=0.5)
				M = 1.0 - Km / N
			return Kunsat = Ks * (Se.^L) * ( 1. - (1. - Se.^ (1.0 / M) ) .^ M ) .^ 2.0
			end # function Se_2_KUNSAT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION vg : ∂K∂Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂K∂Ψ(Ψ₁, iZ, hydroParam, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψvg=hydroParam.Ψvg[iZ], N=hydroParam.N[iZ], Ks=hydroParam.Ks[iZ], Km=hydroParam.Km[iZ], L=0.5)
				M = 1.0 - Km/N
		
				∂K∂Ψ = Ks * (L * (-M * (N * (1 / Ψvg) * (Ψ₁ / Ψvg) ^ (N - 1)) * (1.0 + (Ψ₁ / Ψvg) ^ N) ^ (-M - 1)) * ((1.0 + (Ψ₁ / Ψvg) ^ N) ^ -M) ^ (L - 1)) * (1.0 - (1.0 - ((1.0 + (Ψ₁ / Ψvg) ^ N) ^ -M) ^ (1.0 / M)) ^ M) ^ 2.0 + Ks *
				((1.0 + (Ψ₁ / Ψvg) ^ N) ^ -M) ^ L * (2.0 * -(M * -((1.0 / M) * (-M * (N * (1 / Ψvg) * (Ψ₁ / Ψvg) ^ (N - 1)) * (1.0 + (Ψ₁ / Ψvg) ^ N) ^ (-M - 1)) * ((1.0 + (Ψ₁ / Ψvg) ^ N) ^ -M) ^ (1.0 / M - 1)) * (1.0 - ((1.0 + (Ψ₁ / Ψvg) ^ N) ^ -M) ^ (1.0 / M)) ^ (M - 1)) * (1.0 - (1.0 - ((1.0 + (Ψ₁ / Ψvg) ^ N) ^ -M) ^ (1.0 / M)) ^ M))

				return ∂K∂Ψ
			end # function ∂K∂Ψ

	end #module vg ...............................................


	# =============================================================
	#		MODULE BROOKS AND COOREY
	# =============================================================
	module bc
		import ..option, ..wrc
		export Ψ_2_KUNSAT, Se_2_KUNSAT, ∂K∂Ψ

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION bc : Ψ_2_KUNSAT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function  Ψ_2_KUNSAT(Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψbc=hydroParam.Ψbc[iZ], λbc=hydroParam.λbc[iZ], Ks=hydroParam.Ks[iZ])
				
				M = -3.0 * λbc - 2.0

				if Ψ₁ > Ψbc
					return Kunsat = Ks * (Ψ₁ / Ψbc) ^ M 
				else 
					return Kunsat = Ks
				end

			end #function Ψ_2_KUNSAT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION bc : Se_2_KUNSAT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function  Se_2_KUNSAT(Se, iZ, hydroParam, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψbc=hydroParam.Ψbc[iZ], λbc=hydroParam.λbc[iZ], Ks=hydroParam.Ks[iZ])
				
				M = -3.0 * λbc - 2.0
				return Kunsat = Ks * Se .^ M 
			
			end # function Se_2_KUNSAT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION bc : ∂K∂Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂K∂Ψ(Ψ₁, iZ, hydroParam, θs=hydroParam.θs[iZ], θr=hydroParam.θr[iZ], Ψbc=hydroParam.Ψbc[iZ], λbc=hydroParam.λbc[iZ], Ks=hydroParam.Ks[iZ])
				
				M = -3.0 * λbc - 2.0
		
				return ∂K∂Ψ = Ks * M / Ψbc * (Ψ₁/Ψbc) ^ (M-1.0) 

			end # function ∂K∂Ψ


	end #module bc ...............................................


	# =============================================================
	#		MODULE CLAPP AND HORNBERGER
	# =============================================================
	module ch
		import ..option, ..wrc
		export Ψ_2_KUNSAT, Se_2_KUNSAT, ∂K∂Ψ

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION ch : Ψ_2_KUNSAT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function  Ψ_2_KUNSAT(Ψ₁, iZ, hydroParam; θs=hydroParam.θs[iZ], Ψch=hydroParam.Ψch[iZ], λch=hydroParam.λch[iZ], Ks=hydroParam.Ks[iZ])
				
				M = -3.0 * λch - 2.0

				if Ψ₁ > Ψch
					return Kunsat = Ks * (Ψ₁ / Ψch) ^ M 
				else 
					return Kunsat = Ks
				end

			end #function Ψ_2_KUNSAT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION ch : Se_2_KUNSAT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function  Se_2_KUNSAT(Se, iZ, hydroParam, θs=hydroParam.θs[iZ], Ψch=hydroParam.Ψbc[iZ], λch=hydroParam.λch[iZ], Ks=hydroParam.Ks[iZ])
					
				M = -3.0 * λch - 2.0
				return Kunsat = Ks * Se .^ M 
			
			end # function Se_2_KUNSAT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION ch : ∂K∂Ψ
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∂K∂Ψ(Ψ₁, iZ, hydroParam, θs=hydroParam.θs[iZ], Ψch=hydroParam.Ψch[iZ], λch=hydroParam.λch[iZ], Ks=hydroParam.Ks[iZ])
				
				M = -3.0 * λch - 2.0
		
				return ∂K∂Ψ = Ks * M / Ψch * (Ψ₁/Ψch) ^ (M-1.0) 

			end # function ∂K∂Ψ

	end #module ch ...............................................


end # module kunsat 
# ...........................................................................