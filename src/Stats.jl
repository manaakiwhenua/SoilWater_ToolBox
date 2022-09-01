module stats
	import ..wrc, ..param
	export RMSE, NASH_SUTCLIFE_MINIMIZE, NASH_SUTCLIFFE_θΨ, RELATIVE_ERR, NASH_SUTCLIFFE_EFFICIENCY, LINEAR_REGRESSION
	using Statistics

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : NASH_SUTCLIFE_MINIMIZE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RMSE(Obs, Sim; Power=2.0)
			N = length(Obs)
			Err = 0.0
			iCount = 1
			

			for i = 1:N
				if !(isnan(Sim[i]) || isnan(Obs[i]))
					Err += abs(Sim[i] - Obs[i]) ^ Power
					iCount += 1
				end
			end
		return (Err / Float64(iCount)) ^ (1.0 / Power)
		end # function RMSE


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : NASH_SUTCLIFE_MINIMIZE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function NSE(Obs, Sim; Power=2.0)
			N = length(Obs)
			Err = 0.0
			iCount = 1

			∑Obs = 0.0
			iCount = 1
			for i = 1:N
				if !(isnan(Sim[i]) || isnan(Obs[i]))
					∑Obs += Obs[i]
					iCount += 1
				end
			end
			Obs_Mean = ∑Obs / Float64(iCount)
			
			Obs_Mean_Err = 0.0
			for i = 1:N
				if !(isnan(Sim[i]) || isnan(Obs[i]))
					Err += abs(Sim[i] - Obs[i]) ^ Power

					Obs_Mean_Err += abs(Obs_Mean - Obs[i]) ^ Power
				end
			end
		return 1.0 - (Err / Obs_Mean_Err)
		end # function RMSE



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : NASH_SUTCLIFE_MINIMIZE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function NASH_SUTCLIFE_MINIMIZE(Obs, Sim; Power=2.0)
			N = length(Obs)
			Obs_Mean = Statistics.mean(Obs[1:N])
			Obs_Mean_Err = Statistics.sum(abs.(Obs_Mean .- Obs[1:N]) .^ Power)

			if Obs_Mean_Err < 0.000000000000001   
				Obs_Mean_Err = 1.0
			end
			
			Err = 0.0
			for i = 1:N
				Err += abs(Sim[i] - Obs[i]) ^ Power
			end
		return Err / Obs_Mean_Err 
		end  # function: NASH_SUTCLIFE_MINIMIZE


	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : NASH_SUTCLIFFE_θΨ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function NASH_SUTCLIFFE_θΨ(N_SoilSelect, N_Data, Ψ_Sim, θ_Sim, hydro)
			Nse = zeros(Float64, N_SoilSelect)

			for iZ = 1:N_SoilSelect	
				θΨ = zeros(Float64, N_Data[iZ])
				for iRpart = 1:N_Data[iZ]
					θΨ[iRpart] = wrc.Ψ_2_θDual(Ψ_Sim[iZ,iRpart], iZ, hydro)
				end
				Nse[iZ] = 1.0 - NASH_SUTCLIFE_MINIMIZE(θΨ[1:N_Data[iZ]], θ_Sim[iZ,1:N_Data[iZ]])	
			end

			# Cumulating the objective function to get the overview
			Nse_Mean = round(Statistics.mean(max.(Nse[1:N_SoilSelect],0.0)), digits=3)  # in case of negative value then it is set to 0
			Nse_Std  = round(Statistics.std(max.(Nse[1:N_SoilSelect],0.0)), digits=3)   # in case of negative value then it is set to 0

			return Nse, Nse_Mean, Nse_Std
		end # function NASH_SUTCLIFFE_θΨ

		

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : NASH_SUTCLIFFE_EFFICIENCY
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function NASH_SUTCLIFFE_EFFICIENCY(;Obs=Obs, Sim=Sim, Power=2.0)
			return Nse = 1 - NASH_SUTCLIFE_MINIMIZE(Obs, Sim; Power=Power)
		end  # function: NASH_SUTCLIFFE_EFFICIENCY



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RELATIVE_ERR
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RELATIVE_ERR(; Obs=Obs, Sim=Sim)
			return Err = 1.0 - abs(Obs - Sim) / Obs
		end



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : LINEAR_REGRESSION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function LINEAR_REGRESSION(Xdata, Ydata)
			X = hcat(ones(length(Xdata)),Xdata)
			Y = Ydata
			return Intercept, Slope = inv(X'*X)*(X'*Y)
		end # function LINEAR_REGRESSION		function RELATIVE_ERR(;Obs=Obs, Sim=Sim)

end # module