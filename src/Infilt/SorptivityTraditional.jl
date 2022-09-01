# =============================================================
#		MODULE: sorptivity
# =============================================================
module sorptivity
	import ..wrc, ..kunsat, ..option
	import QuadGK
	export SORPTIVITY, DIFFUSIVITY

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : SORPTIVITY
			# https://juliamath.github.io/QuadGK.jl/latest/
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function SORPTIVITY(θ_Ini, iZ, hydro; Rtol=10^-1.0) #-4 Rtol=10^-7.0

			function SORPTIVITY_FUNC(θ, θ_Ini, iZ, hydro)

				Diffusivity = DIFFUSIVITY(θ, iZ, hydro)

				if option.infilt.SorptivityModel == :Parlange # <>=<>=<>=<>=<> strong non-linear behaviur for diffusivity 
					return Sorptivity = (hydro.θs[iZ] + θ - 2.0 * θ_Ini) * Diffusivity

				elseif  option.infilt.SorptivityModel == :Brutsaert  # <>=<>=<>=<>=<>  maximum
					return Sorptivity = 2.0 * ( ((hydro.θs[iZ] -  θ_Ini)^ 0.5) * (θ - θ_Ini)^ 0.5   ) * Diffusivity

				elseif  option.infilt.SorptivityModel == "Option3"  # <>=<>=<>=<>=<> Flux concentration function = 1
					return Sorptivity = 2.0 * (  (θ - θ_Ini)   ) * Diffusivity

				elseif  option.infilt.SorptivityModel == "Simmons&McKeon"  # Worst <>=<>=<>=<>=<> intermediate
					return Sorptivity = 2.0 * (  (θ - θ_Ini) / ((θ - θ_Ini) / (hydro.θs[iZ] - θ_Ini)) ^ (2.0 - 4.0 * π) ) * Diffusivity

				elseif  option.infilt.SorptivityModel == :Philip_Knight  # Best <>=<>=<>=<>=<>  minimum
					return Sorptivity = 2.0 * (  (hydro.θs[iZ] - θ_Ini)  ) * Diffusivity

				elseif  option.infilt.SorptivityModel == "Whiteetal"  #  <>=<>=<>=<>=<>  maximum
					return Sorptivity = 2.0 * (  (θ  - θ_Ini) / sin((π/2)*((θ - θ_Ini) / (hydro.θs[iZ] - θ_Ini))^(π/4) ) ) * Diffusivity
				end # option.infilt
			end # function: SORPTIVITY_FUNC

			return ( QuadGK.quadgk(θ -> SORPTIVITY_FUNC(θ, θ_Ini, iZ, hydro), θ_Ini, hydro.θs[iZ] - eps(), rtol=Rtol )[1]) ^ 0.5  
		end  # function: SORPTIVITY_MODEL


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : DIFFUSIVITY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function DIFFUSIVITY(θ, iZ, hydro; Diffusivity_Min=10^-5.0) #-5		
			Kunsat = kunsat.θ_2_KUNSAT(θ, iZ, hydro)
			Ψ = wrc.θ_2_ΨDual(θ, iZ, hydro)
			∂θ∂Ψ = wrc.∂θ∂Ψ(Ψ, iZ, hydro)

			if ∂θ∂Ψ >  Diffusivity_Min
				return Diffusivity = Kunsat / ∂θ∂Ψ
			else
				return Diffusivity = 0.0
			end
		end  # function: DIFFUSIVITY

	
end  # module: sorptivity
# ............................................................