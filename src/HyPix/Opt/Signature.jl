# =============================================================
#		module: Signature
# =============================================================
module signature

	import ..param, ..tool, ..wrc
	import Dates: value, DateTime, year, month, day, hour, minute, second

	export SIGNATURE

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : SIGNATURE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	function SIGNATURE(∑T, obsθ, hydroHorizon, N_iRoot, N_iT, N_iZ, veg, ΔRootDensity, Ψ)

		# SIMULATED 
			# Simulated which matches time of θobs
				N_∑T_Plot, True = tool.readWrite.DATA_2_ΔTnew(∑T, N_iT, obsθ.∑T[1:obsθ.N_iT]) # time resolution of θsim at time resolution obsθ
	
				Ψ_Sim = Ψ[True[1:N_iT],1:N_iZ]

			# Simulated which matches the number of layers as θobs
            ΔT_Plot   = Matrix{Float64}(undef, N_∑T_Plot, N_iZ) # Reserving memory
            Ψ_SimReduced = Matrix{Float64}(undef, N_∑T_Plot, N_iZ)

				ΔT_Plot[1] = obsθ.∑T[1]
				for iT = 1:N_∑T_Plot
					for iDepth = 1:obsθ.Ndepth
						Ψ_SimReduced[iT,iDepth] = Ψ_Sim[iT, obsθ.ithetaObs[iDepth]]
					end  # for obsθ.Ndepth
					
					if iT ≥ 2
						ΔT_Plot[iT] = obsθ.∑T[iT] - obsθ.∑T[iT-1]
					end
				end  # for iT = 1:N_∑T_Plot

			
		# OBSERVED
			# New ΔRootDensity which has the number f layers as observed θ
				ΔRootDensity_Plot = zeros(Float64, obsθ.Ndepth) # Reserving memory

				iDepth = 1
				for iZ=1:N_iRoot
					if obsθ.ithetaObs[iDepth] < iZ	
						if iDepth ≤ obsθ.Ndepth
							iDepth += 1
						end
					end #  obsθ.ithetaObs[iDepth] ≤ iZ
					ΔRootDensity_Plot[iDepth] += ΔRootDensity[iZ]
				end # for iZ=1:N_iRoot
			
			# Ψ_Obs	
				Ψ_Obs = Matrix{Float64}(undef, N_∑T_Plot, obsθ.Ndepth) # Reserving memory

				for iT=1:N_∑T_Plot, iZ=1:obsθ.Ndepth
					Ψ_Obs[iT,iZ] = wrc.θ_2_ΨDual(min.( max.(obsθ.θobs[iT,iZ], hydroHorizon.θr[iZ]), hydroHorizon.θs[iZ]), iZ, hydroHorizon)
				end # for iT iZ

				Signature_Deficit_Obs, Signature_Max_Obs, Signature_Saturated_Obs, Signature_Senescence_Obs, Signature_Deficit_Sim, Signature_Max_Sim, Signature_Saturated_Sim, Signature_Senescence_Sim = SIGNATURE_Ψ(obsθ, N_∑T_Plot, veg, ΔRootDensity_Plot, ΔT_Plot, Ψ_Obs, Ψ_SimReduced)

		end # function Signatures

	#= ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			FUNCTION : Signature_Ψ
	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ =#
		function SIGNATURE_Ψ(obsθ, N_∑T_Plot, veg, ΔRootDensity_Plot, ΔT_Plot, Ψ_Obs, Ψ_SimReduced; Wsignature=0.75)
       
			# OBSERVED
			Signature_Deficit_Obs, Signature_Max_Obs, Signature_Saturated_Obs, Signature_Senescence_Obs, ∑ΔT_Month = SIGNATURE_FEDDES(obsθ, N_∑T_Plot, veg, ΔRootDensity_Plot, ΔT_Plot, Ψ_Obs)
			
			# SIMULATED
			Signature_Deficit_Sim, Signature_Max_Sim, Signature_Saturated_Sim, Signature_Senescence_Sim, ~ = SIGNATURE_FEDDES(obsθ, N_∑T_Plot, veg, ΔRootDensity_Plot, ΔT_Plot, Ψ_SimReduced)
		
			return Signature_Deficit_Obs, Signature_Max_Obs, Signature_Saturated_Obs, Signature_Senescence_Obs, Signature_Deficit_Sim, Signature_Max_Sim, Signature_Saturated_Sim, Signature_Senescence_Sim

		end  # function: SIGNATURE_Ψ


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : SIGNATURE_FEDDES
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function SIGNATURE_FEDDES(obsθ, N_∑T_Plot, veg, ΔRootDensity_Plot, ΔT_Plot, Ψ₂)

         Signature_Deficit    = zeros(Float64, 12)
         Signature_Max        = zeros(Float64, 12)
         Signature_Saturated  = zeros(Float64, 12)
         Signature_Senescence = zeros(Float64, 12)
         ∑ΔT_Month            = zeros(Float64, 12)
		
			for iT=2:N_∑T_Plot
				iMonth = Int64(month(obsθ.Date[iT]))

            Signature_Deficit₁    = 0.0
            Signature_Max₁        = 0.0
            Signature_Saturated₁  = 0.0
            Signature_Senescence₁ = 0.0

				for iZ=1:obsθ.Ndepth
					# Signature_Saturated
						if Ψ₂[iT,iZ] < veg.Ψfeddes2
							Signature_Saturated₁ += ΔRootDensity_Plot[iZ] # Signature_Saturated = Signature_Saturated + ΔRootDensity_Plot[iZ]
											
					# Signature_Max
						elseif veg.Ψfeddes2 ≤ Ψ₂[iT,iZ] ≤ veg.Ψfeddes3
							Signature_Max₁ += ΔRootDensity_Plot[iZ]
						
					# Signature_Deficit
						elseif veg.Ψfeddes3 < Ψ₂[iT,iZ] < veg.Ψfeddes4 
							Signature_Deficit₁ += ΔRootDensity_Plot[iZ]

					# Signature_Senescence
						elseif veg.Ψfeddes4 ≤ Ψ₂[iT,iZ]
							Signature_Senescence₁ += ΔRootDensity_Plot[iZ]
						end # if
				end # for iZ=1:N_iZ₂

				# Giving a UNIT OF time
					Signature_Deficit[iMonth] 		= Signature_Deficit₁ * ΔT_Plot[iT]
					Signature_Max[iMonth] 			= Signature_Max₁ * ΔT_Plot[iT]
					Signature_Saturated[iMonth] 	= Signature_Saturated₁ * ΔT_Plot[iT]
					Signature_Senescence[iMonth] 	= Signature_Senescence₁ * ΔT_Plot[iT]

				# Time spend every month
					∑ΔT_Month[iMonth] += ΔT_Plot[iT] # giving the time spend at every month

			end # for iT=1:N_∑T_Plot

			# NORMALISING
				for iMonth=1:12
					if ∑ΔT_Month[iMonth] > 1
                  Signature_Deficit[iMonth]		/= ∑ΔT_Month[iMonth]
                  Signature_Max[iMonth] 			/= ∑ΔT_Month[iMonth]
                  Signature_Saturated[iMonth] 	/= ∑ΔT_Month[iMonth]
                  Signature_Senescence[iMonth] 	/= ∑ΔT_Month[iMonth]
					else
                  Signature_Deficit[iMonth]    = -999
                  Signature_Max[iMonth]        = -999
                  Signature_Saturated[iMonth]  = -999
                  Signature_Senescence[iMonth] = -999
					end
				end

			return Signature_Deficit, Signature_Max, Signature_Saturated, Signature_Senescence, ∑ΔT_Month

		end  # function: SIGNATURE_FEDDES


end  # module: Signature
# ............................................................