# =============================================================
#		MODULE: sorptivity
# =============================================================
module sorptivity
	import ..wrc, ..kunsat, ..option
	import QuadGK
	import SpecialFunctions: erfc, erfcinv
	export SORPTIVITY

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : SORPTIVITY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function SORPTIVITY(θ_Ini, iZ, hydroInfilt; Rtol= 10^-8.0, SorptivityModelScaled = true) #Rtol=10^-3.0

			# INITIALIZING
				Ψ_Sat = 0.0

				Se_Ini = wrc.θ_2_Se(θ_Ini, iZ, hydroInfilt)

				Ψ_Ini = wrc.θ_2_ΨDual(θ_Ini, iZ, hydroInfilt)

				SeIni_⬙ = (1.0 + Se_Ini) / 2.0

				θ⬙ = wrc.Se_2_θ(SeIni_⬙, iZ, hydroInfilt)

				Ψ⬙ = wrc.θ_2_ΨDual(θ⬙, iZ, hydroInfilt)

		
			if option.infilt.SorptivitySplitModel == :Split  # <>=<>=<>=<>=<>
				# Sorptivity based on θ
					function SORPTIVITY_θ²(hydroInfilt, iZ, θ, θ_Ini)
						Se = wrc.θ_2_Se(θ, iZ, hydroInfilt)

						return Sorptivity_θ = DIFFUSIVITY_θ(θ, iZ, hydroInfilt) * (hydroInfilt.θs[iZ] + θ - 2.0 * θ_Ini)
					end # SORPTIVITY_θ² ~~~~~~~~~~~~~~~~~

					Sorptivity_θ² = QuadGK.quadgk(θ -> SORPTIVITY_θ²(hydroInfilt, iZ, θ, θ_Ini), θ_Ini, θ⬙, rtol=Rtol)[1]

				# Sorptivity based on Ψ
					function SORPTIVITY_Ψ²(hydroInfilt, iZ, θ_Ini, Ψ)
						θ = wrc.Ψ_2_θDual(Ψ, iZ, hydroInfilt)
						return Sorptivity_Ψ = kunsat.Ψ_2_KUNSAT(Ψ, iZ, hydroInfilt) * (hydroInfilt.θs[iZ] + θ - 2.0 * θ_Ini)
					end # SORPTIVITY_Ψ² ~~~~~~~~~~~~~~~~~

					Sorptivity_Ψ² = QuadGK.quadgk(Ψ -> SORPTIVITY_Ψ²(hydroInfilt, iZ, θ_Ini, Ψ), Ψ_Sat, Ψ⬙, rtol=Rtol)[1]

			return Sorptivity = √(max(Sorptivity_θ², eps()) + max(Sorptivity_Ψ², eps()))
			elseif option.infilt.SorptivitySplitModel == :Split_η # <>=<>=<>=<>=<>
	
				Ψ⬙_η = wrc.θ_2_ΨDual(θ⬙, iZ, hydroInfilt) / hydroInfilt.Ψm[iZ] # dimensionless water potential

				ΨSat_η = Ψ_Sat / hydroInfilt.Ψm[iZ]

				# Scaled sorptivity based on Se
					function SORPTIVITY_Se_η(Se, hydroInfilt, iZ, Se_Ini)
						return Sorptivity_Se_η = DIFFUSIVITY_Se_η(Se, iZ, hydroInfilt) * FLUXCONC(Se, Se_Ini)
					end # SORPTIVITY_Se_η

					Sorptivity_Se_η = QuadGK.quadgk(Se -> SORPTIVITY_Se_η(Se, hydroInfilt, iZ, Se_Ini), Se_Ini, SeIni_⬙, rtol=Rtol)[1]
				
				# Scaled SORPTIVITY_Se_Ψ
					function SORPTIVITY_Se_Ψ_η(Ψ_η, hydroInfilt, iZ, Se_Ini)
						Ψ₂ = Ψ_η * hydroInfilt.Ψm[iZ]

						Se = wrc.Ψ_2_SeDual(Ψ₂, iZ, hydroInfilt)

						return Sorptivity_Se_Ψ_η = kunsat.Se_2_KR(Se, iZ, hydroInfilt) * FLUXCONC(Se, Se_Ini)
					end # SORPTIVITY_Se_Ψ

					Sorptivity_Se_Ψ_η = QuadGK.quadgk(Ψ_η -> SORPTIVITY_Se_Ψ_η(Ψ_η, hydroInfilt, iZ, Se_Ini), ΨSat_η, Ψ⬙_η, rtol=Rtol)[1]
	
				# To avoid numerical instability when Se = 1
				Sorptivity_η = √(max(Sorptivity_Se_η , eps()) + max(Sorptivity_Se_Ψ_η , eps()))
	
			return Sorptivity_Scaled = Sorptivity_η * √(hydroInfilt.Ks[iZ] * (hydroInfilt.θs[iZ] - hydroInfilt.θr[iZ]) * hydroInfilt.Ψm[iZ])
			end # option.infilt.SorptivitySplitModel

		end # function SORPTIVITY


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#  	FUNCTION : DIFFUSIVITY
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function DIFFUSIVITY_θ(θ, iZ, hydroInfilt)
				Kunsat = kunsat.θ_2_KUNSAT(θ, iZ, hydroInfilt)

				Ψ = wrc.θ_2_ΨDual(θ, iZ, hydroInfilt)
				
				∂Ψ∂θ = -wrc.∂Ψ∂θ(θ, iZ, hydroInfilt)
				
				return Diffusivity = Kunsat * ∂Ψ∂θ
			end  # function: DIFFUSIVITY_θ ~~~~~~~~~~~~~~~~~


			# DIFFUSIVITY_Se_η
			function DIFFUSIVITY_Se_η(Se, iZ, hydroInfilt)
				Kr = kunsat.Se_2_KR(Se, iZ, hydroInfilt)

				Ψ = wrc.Se_2_ΨDual(Se, iZ, hydroInfilt)

				∂Ψ∂Se_η = -wrc.∂Ψ∂Se(Se, iZ, hydroInfilt) / hydroInfilt.Ψm[iZ]# negative sign needed because Ψ is set positive

				return Diffusivity = Kr * ∂Ψ∂Se_η
			end  # function: DIFFUSIVITY_Se_η ~~~~~~~~~~~~~~~~~


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : FLUXXCONCENTRATION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function FLUXCONC(Se, Se_Ini)
			if option.infilt.SorptivityModel == :Parlange # <>=<>=<>=<>=<> soils with strong non-linear diffusivity behaviors
				return FluxConc = (1.0 + Se - 2.0 * Se_Ini) # to avoid numerical problems it is set analyticaly

			elseif option.infilt.SorptivityModel == :Crank # <>=<>=<>=<>=<> linear soil with constant diffusivity
				return FluxConc = (2.0 * (Se - Se_Ini)) / exp(-(erfcinv(max(Se - Se_Ini, eps()) / (1.0 - Se_Ini)))^2 )

			elseif option.infilt.SorptivityModel == :Philip_Knight # <>=<>=<>=<>=<> soils with Dirac δ-function diffusivity (i.e. Green and Ampt model)
				return FluxConc = (2.0 * (1.0 - Se_Ini)) # to avoid numerical problems it is set analyticaly 

			elseif option.infilt.SorptivityModel == :Brutsaert # <>=<>=<>=<>=<> soil with moderate non-linear diffusivity behaviors
				return FluxConc = (2.0 * √((Se - Se_Ini) * (1.0 - Se_Ini)))  # to avoid numerical problems it is set analyticaly
			end # option.infilt.SorptivityModel
		end  # function: FLUXXCONCENTRATION

end  # module: sorptivity
# ............................................................