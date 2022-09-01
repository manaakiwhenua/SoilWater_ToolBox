module climate
	import Dates: value, DateTime
	import ..option

	export CLIMATE

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :   CLIMATE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function CLIMATE(clim)

			fill(0.0::Float64, clim.N_Climate)

         ∑Pr_Climate      = fill(0.0::Float64, clim.N_Climate)
         ∑Pet_Climate	  = fill(0.0::Float64, clim.N_Climate)
         ∑T_Climate       = fill(0.0::Float64, clim.N_Climate)
         Temp             = fill(0.0::Float64, clim.N_Climate)

			# Computing Pr, PotEvap, Temp in the correct format
	 		Start_Date = clim.Date[1]

			 # Taking into acount that ΔT is the difference of time of T[iT]-T[iT-1]
          ∑Pr_Climate[1]  = 0.0
          ∑Pet_Climate[1] = 0.0
          ∑T_Climate[1]   = 0.0
			 Temp[1]         = 0.0
			 
			 for iT = 2:clim.N_Climate
				#Computing cumulative time 
				∑T_Climate[iT] = value(clim.Date[iT] - Start_Date) / 1000

				# Cumulative
				if !(option.hyPix.RainfallInterception)
					∑Pr_Climate[iT] = ∑Pr_Climate[iT-1] + clim.Pr[iT]  # Cumulate Pr
					∑Pet_Climate[iT] = ∑Pet_Climate[iT-1] + clim.Pet[iT] # Cumulative Potential evaporation
				end
				
				# Tempoerature no change 
				Temp[iT] = clim.Temp[iT]
			end # for

			N_∑T_Climate = Int(∑T_Climate[clim.N_Climate])

		return ∑Pet_Climate, ∑Pr_Climate, ∑T_Climate, N_∑T_Climate, Temp
		end # function CLIMATE

end # module climate
# ...........................................................................................