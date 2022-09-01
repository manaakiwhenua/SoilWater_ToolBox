# =============================================================
#		MODULE: waterbalance
# =============================================================
module waterBalance

	export WATERBALANCE

	function WATERBALANCE(∑T, obsθ, discret, hydro, N_iRoot::Int64, N_iT::Int64, N_iZ::Int64, Q, ΔSink, ΔT, θ, Ψ)
      ∑WaterBalance_η = fill(0.0::Float64, N_iT)
      ∑ΔSink          = fill(0.0::Float64, N_iT)
      ∑ΔQtop          = 0.0 ::Float64
      ∑ΔQbottom       = 0.0 ::Float64
      ∑∑WaterBalance  = 0.0 ::Float64
		ΔStorage        = 0.0 ::Float64
		ΔStorageSo =  0.0 ::Float64
      i∑T_CalibrStart    = 1::Int64

		# Starting to compute the waterbalance after the warmup period
			while ∑T[i∑T_CalibrStart] < obsθ.∑T[1]
				i∑T_CalibrStart += 1
			end
		
		for iT=i∑T_CalibrStart:N_iT
			# Computing ΔStorage
				ΔStorage = 0.0

				@fastmath @inbounds for iZ = 1:N_iZ
					ΔStorage += discret.ΔZ[iZ] * ( (θ[iT,iZ] - θ[i∑T_CalibrStart-1,iZ]) )

					# ΔStorage += discret.ΔZ[iZ] * ( (θ[iT,iZ] - θ[i∑T_CalibrStart-1,iZ]) - hydro.So[iZ] * (Ψ[iT,iZ] - Ψ[iT-1,iZ]) * (θ[iT,iZ] / hydro.θs[iZ]) )
					
					ΔStorageSo += discret.ΔZ[iZ] * hydro.So[iZ] * (Ψ[iT,iZ] - Ψ[iT-1,iZ]) * (θ[iT,iZ] / hydro.θs[iZ])
				end # for iT=1:N_iZ

			# Sink term
				∑ΔSink[iT] = ∑ΔSink[iT-1]
				@fastmath @inbounds @simd for iZ = 1:N_iRoot
					∑ΔSink[iT] += ΔSink[iT,iZ]
				end # for: iZ = 1:N_iZ

			# Cumulative water entering top cell
				∑ΔQtop += ΔT[iT] * Q[iT,1]

			# Cumulative water leaving bottom cell
				∑ΔQbottom += ΔT[iT] * Q[iT,N_iZ+1]

			∑∑WaterBalance = ΔStorage - (∑ΔQtop - ∑ΔQbottom) + ∑ΔSink[iT] - ΔStorageSo

			∑WaterBalance_η[iT] = ∑∑WaterBalance / ∑ΔQtop
		end  # for iT=1:N_iT

		println("	=== ΔStorage = $ΔStorage")

		return ∑∑WaterBalance, ∑WaterBalance_η, ∑ΔSink, i∑T_CalibrStart, ΔStorage
	end  # function: WATERBALANCE

end  # module waterbalance
# ............................................................