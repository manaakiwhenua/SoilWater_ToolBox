# =============================================================
#		module: discret
# =============================================================
module discretization
	import ..param
	export DISCRETIZATION

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :  DISCRETIZATION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		struct DISCRET
			N_iZ
			Z_CellUp
			Znode
			ΔZ
			ΔZ_⬓
			ΔZ_Aver
			ΔZ_W
		end
		
		function DISCRETIZATION(N_iZ::Int64, Z)
			Z_CellUp = fill(0.0::Float64, N_iZ)
			Znode    = fill(0.0::Float64, N_iZ)
			ΔZ       = fill(0.0::Float64, N_iZ)
			ΔZ_⬓     = fill(0.0::Float64, N_iZ)
			ΔZ_Aver  = fill(0.0::Float64, N_iZ+1)
			ΔZ_W     = fill(0.0::Float64, N_iZ+1)
			
			# Cell 1
				ΔZ[1]       = Z[1]
				ΔZ_⬓[1]     = ΔZ[1] * 0.5
				Z_CellUp[1] = 0.0
				Znode[1]    = ΔZ_⬓[1]
				ΔZ_Aver[1]  = ΔZ_⬓[1]
				ΔZ_W[1]     = 1.0

			# All Cells
			for iZ = 2:N_iZ
				ΔZ[iZ]       = Z[iZ] - Z[iZ-1]

				ΔZ_⬓[iZ]    = ΔZ[iZ] * 0.5

				Znode[iZ]    = Z[iZ] - ΔZ_⬓[iZ]

				ΔZ_Aver[iZ]  = (ΔZ[iZ] + ΔZ[iZ-1]) * 0.5

				ΔZ_W[iZ]     = ΔZ[iZ] /  (ΔZ[iZ] + ΔZ[iZ-1])

				Z_CellUp[iZ] = Z[iZ] - ΔZ[iZ]
			end # for

			ΔZ_Aver[N_iZ+1] = ΔZ[N_iZ] * 0.5

			ΔZ_W[N_iZ+1] = 0.0

		return discret = DISCRET(N_iZ, Z_CellUp, Znode, ΔZ, ΔZ_⬓, ΔZ_Aver, ΔZ_W)
		end # function DISCRETIZATION


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : DISCRETISATION_AUTO
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	""" 
					DISCRETISATION_AUTO(Nlayer, Zlayer)

	Automatically performs the discretisatio of the HyPix model wheh you enter the depth of the layers
	"""
		function DISCRETISATION_AUTO(; Nlayer, Zlayer, Zroot, θᵢₙᵢ)

			ΔZlayer = fill(0.0::Float64, Nlayer)

			# Computing ΔZlayer
				ΔZlayer[1]= Zlayer[1]
				for iZ = 2:Nlayer
					ΔZlayer[iZ] = Zlayer[iZ] - Zlayer[iZ-1]
				end # for

			# Computing the number of discretization
            ΔZcell    = []
            Layer     = []
            θᵢₙᵢ_Cell = []
				for iLayer =1:Nlayer

					if  Zlayer[iLayer] < Zroot
						ΔZ_Max = param.hyPix.ΔZrz_Max
					else
						ΔZ_Max = param.hyPix.ΔZdeep_max
					end

					Nsplit = ceil(ΔZlayer[iLayer] / ΔZ_Max) # Number of splitting from Layer->Cell
					ΔZcell₀ = ΔZlayer[iLayer] / Float64(Nsplit)

					for iDiscret=1:Nsplit
						append!(ΔZcell, ΔZcell₀)
                  append!(Layer, iLayer)
						append!(θᵢₙᵢ_Cell, θᵢₙᵢ[iLayer])
					end
				end # ilayer
				N = length(ΔZcell)

			# Computing the ∑ΔZcell
				Z = fill(0.0::Float64, N)

				Z[1] = ΔZcell[1]
				for iZ=2:N
					Z[iZ] = Z[iZ-1] + ΔZcell[iZ]
				end
			
		return Layer, Z, θᵢₙᵢ_Cell
		end  # function: DISCRETISATION_AUTO
	
end  # module: discret
# ............................................................