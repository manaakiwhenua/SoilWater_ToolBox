module OfHydrolab 
	import ..option, ..stats, ..wrc, ..kunsat
	export  OF_WRC_KUNSAT

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OF_WRC_KUNSAT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
		function OF_WRC_KUNSAT(iZ, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro, optim, optionHydro; Wof = 0.5) 

			# === OF θΨ ====
				θ_Obs = Array{Float64}(undef, N_θΨ[iZ])
				θ_Sim = Array{Float64}(undef, N_θΨ[iZ])

				for iΨ = 1:N_θΨ[iZ]
					θ_Obs[iΨ] = θ_θΨ[iZ,iΨ]
					θ_Sim[iΨ] = wrc.Ψ_2_θDual(Ψ_θΨ[iZ,iΨ], iZ, hydro)
				end # for iΨ = 1:N_θΨ[iZ]

				Of_θΨ = stats.NASH_SUTCLIFE_MINIMIZE(θ_Obs[1:N_θΨ[iZ]], θ_Sim[1:N_θΨ[iZ]])

			# === OF Kunsat ====
			if optionHydro.KunsatΨ || optionHydro.Kunsat_JustRun
				if "Ks" ∈ optim.ParamOpt
					iStart = 1
				else
					iStart = 2
				end

				# Of_Kunsat = 0.0
				Kunsat_Obs_Ln = Array{Float64}(undef, N_KΨ[iZ])
				Kunsat_Sim_Ln = Array{Float64}(undef, N_KΨ[iZ])

				for iΨ = iStart:N_KΨ[iZ]
					Kunsat_Obs_Ln[iΨ] = log1p(K_KΨ[iZ,iΨ])
						
					Kunsat_Sim_Ln[iΨ] = log1p(kunsat.Ψ_2_KUNSAT(Ψ_KΨ[iZ,iΨ], iZ, hydro))
				end # for iΨ = 1:N_KΨ[iZ]

				Of_Kunsat = stats.NASH_SUTCLIFE_MINIMIZE(Kunsat_Obs_Ln[iStart:N_KΨ[iZ]], Kunsat_Sim_Ln[iStart:N_KΨ[iZ]])			

				if option.hydro.Kunsat_JustRun
					Of = Of_θΨ
				else
					Of = Wof * Of_θΨ + (1.0 - Wof) * Of_Kunsat
				end
			else		
				Of = Of_θΨ
				Of_Kunsat = 0.0
			end #  optionHydro.KunsatΨ

		return Of, Of_θΨ, Of_Kunsat
		end # function OF_WRC_KUNSAT



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : OF_RMSE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
		function OF_RMSE(iZ, θ_θΨ, Ψ_θΨ, N_θΨ, K_KΨ, Ψ_KΨ, N_KΨ, hydro, optim, optionHydro) 

		# === OF θΨ ====
			θ_Obs = Array{Float64}(undef, N_θΨ[iZ])
			θ_Sim = Array{Float64}(undef, N_θΨ[iZ])

			for iΨ = 1:N_θΨ[iZ]
				θ_Obs[iΨ] = θ_θΨ[iZ,iΨ]
				θ_Sim[iΨ] = wrc.Ψ_2_θDual(Ψ_θΨ[iZ,iΨ], iZ, hydro)
			end # for iΨ = 1:N_θΨ[iZ]

			Rmse_θΨ = stats.RMSE(θ_Obs[1:N_θΨ[iZ]], θ_Sim[1:N_θΨ[iZ]])

		# === OF Kunsat ====
			if optionHydro.KunsatΨ ||optionHydro.Kunsat_JustRun
				if  "Ks" ∈ optim.ParamOpt
					iStart = 1
				else
					iStart = 2
				end

				Kunsat_Obs_Ln = Array{Float64}(undef, N_KΨ[iZ])
				Kunsat_Sim_Ln = Array{Float64}(undef, N_KΨ[iZ])

				for iΨ = iStart:N_KΨ[iZ]
					Kunsat_Obs_Ln[iΨ] = log1p(K_KΨ[iZ,iΨ])
						
					Kunsat_Sim_Ln[iΨ] = log1p(kunsat.Ψ_2_KUNSAT(Ψ_KΨ[iZ,iΨ], iZ, hydro))
				end # for iΨ = 1:N_KΨ[iZ]

				Rmse_KΨ = stats.RMSE(Kunsat_Obs_Ln[iStart:N_KΨ[iZ]], Kunsat_Sim_Ln[iStart:N_KΨ[iZ]])

				Rmse = (Rmse_θΨ + Rmse_KΨ) * 0.5
			else		
				Rmse = Rmse_θΨ
				Rmse_KΨ = 0.0
			end #  optionHydro.KunsatΨ

	return Rmse, Rmse_KΨ, Rmse_θΨ
	end # OF_RMSE

end # module ofHydaulic
# ............................................................