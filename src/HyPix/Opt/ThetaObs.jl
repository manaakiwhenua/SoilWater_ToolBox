# =============================================================
#		module: thetaObs
# =============================================================
module thetaObs
	import Dates: value, DateTime
	export thetaObs
	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ΘOBS
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ΘOBS(obsθ, clim, discret, Z)

			# CHECKING DATA CONSISTENCY
				if obsθ.Date[1] < clim.Date[2]
					error("\n Hypix error: Starting date of obsθ  $(obsθ.Date[1]) < starting date of climate data $(clim.Date[1])")
				end # Error checking

				# Checking the celLs
				if obsθ.Z[obsθ.Ndepth] > Z[discret.N_iZ]
					error("\n Hypix error: depth of measured θ deeper than the max depth of discretisation: obsθ.Z[obsθ.Ndepth] > discret.Z[discret.N_iZ]") 
				end

			# COMPUTING CUMULATIVE TIME
				Start_Date = clim.Date[1] 
				for iT=1:obsθ.N_iT
					obsθ.∑T[iT] = value(obsθ.Date[iT] - Start_Date) / 1000
				end  # for it=1:obsθ.N_iT

				# TRANSFORM THE DEPTH OF MEASURED Θ -> CELL DEPTH
				for iDepth = 1:obsθ.Ndepth
					for iZ = 1:discret.N_iZ
						if iZ == 1
							if 0.0 ≤ obsθ.Z[iDepth] ≤ Z[1]
								obsθ.ithetaObs[iDepth] = 1
								break  
							end  # if: discret.Z_CellUp
						elseif iZ ≠ 1
							if Z[iZ-1] ≤ obsθ.Z[iDepth] ≤ Z[iZ]
								obsθ.ithetaObs[iDepth] = iZ
								break  
							end  # if: discret.Z_CellUp
						end # if iZ == 1
					end # iZ = 2:discret.N_iZ						
				end  # for iDepth = 1:obsθ.Ndepth
		
		return obsθ
		end  # function: θOBS

end # module: thetaObs