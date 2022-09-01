# =============================================================
#		MODULE: table
# =============================================================
module table
	import ..path
	import Tables, CSV
	export TABLE_ID

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : TABLE_EXTRAPOINTS_K
	# 		Tabular values of the PSD model
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function TABLE_ID(N_SoilSelect::Int64, Path::String)
			println("    ~  $(Path) ~")

			Id_Select = collect(1:1:N_SoilSelect)

			Select = fill(1::Int64, N_SoilSelect)

			FieldName_String = ["Id", path.Select]

			Output = Tables.table( [Id_Select[1:N_SoilSelect] Select[1:N_SoilSelect]] )
			
			CSV.write(Path, Output, header=FieldName_String, delim=',')
			
		return nothing
		end  # function:  TABLE_ID

	# =============================================================
	#		MODULE: hydroLab
	# =============================================================
	module hydroLab
		import ...path, ...tool, ...param, ...wrc, ...kunsat
		import DelimitedFiles, Tables, CSV
		export θΨK

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : θΨK
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function θΨK(hydro, hydroOther, Id_Select, KunsatModel_Lab, N_SoilSelect)
				println("    ~  $(path.Table_θΨK) ~")

				Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(N_SoilSelect, hydro)

				Matrix2, FieldName_String2 = tool.readWrite.STRUCT_2_FIELDNAME(N_SoilSelect, hydroOther)

				# Concatenating matrices
				Matrix = hcat(Matrix, Matrix2)

				Matrix = hcat(Matrix, KunsatModel_Lab)

				FieldName_String = vcat("Id", FieldName_String, FieldName_String2, "KunsatModel")

				Output = Tables.table([Id_Select Matrix])

				CSV.write(path.Table_θΨK, Output, header=FieldName_String, delim=',' )
			return nothing
			end  # function:  θΨK


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : TABLE_θΨK
		# 		Tabular values of the PSD model
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function TABLE_EXTRAPOINTS_θΨ(hydro, Id_Select, N_SoilSelect, Path, Ψ_Table; Orientation="Horizontal")
				println("    ~  $(Path) ~")

				N_Ψ = Int64(length(Ψ_Table))

				if Orientation == "Horizontal" # <>=<>=<>=<>=<>=<>
					# Writting the Header
						FieldName_String = Vector{String}(undef, (N_Ψ))
						for i =1:N_Ψ
							FieldName_String[i] = string(Int64(Ψ_Table[i]) ) * "mm"
						end
						pushfirst!(FieldName_String, string("Id")) # Write the "Id" at the very begenning
							
				# Computing θ at required θ
					θ₂ = fill(0.0::Float64, (N_SoilSelect, N_Ψ))


					for iZ=1:N_SoilSelect
						for iΨ =1:N_Ψ
							Ψ₂ = Ψ_Table[iΨ]
							θ₂[iZ, iΨ] = wrc.Ψ_2_θDual(Ψ₂, iZ, hydro)
						end # iΨ
					end # iZ

					CSV.write(Path, Tables.table([string.(Id_Select) θ₂]), header=FieldName_String )

				elseif Orientation == "Vertical" # <>=<>=<>=<>=<>=<>
					FieldName_String = ["Id","H[mm]","Theta"]
					N = N_Ψ * N_SoilSelect
               Id₂ = Vector{Int64}(undef, N)
               Ψ₂  = Vector{Float64}(undef,  N)
               θ₂  = Vector{Float64}(undef, N)
					iCount = 1

					for iZ=1:N_SoilSelect
						for iΨ =1:N_Ψ
							Id₂[iCount] = Id_Select[iZ]
							Ψ₂[iCount] = Ψ_Table[iΨ]
							θ₂[iCount] = wrc.Ψ_2_θDual(Ψ₂[iCount], iZ, hydro)
	
							iCount+=1
						end # iΨ
					end # iZ

					Output = Tables.table([string.(Id₂[1:N]) Ψ₂[1:N] θ₂[1:N]])
					CSV.write(Path, Output, header=FieldName_String, delim=',')

				else
					error("SoilWaterToolBox Error in TABLE_EXTRAPOINTS_θΨ $Orientation not understood")
				end

		return nothing
		end  # function:  θΨ


		
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : TABLE_EXTRAPOINTS_K
		# 		Tabular values of the PSD model
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function TABLE_EXTRAPOINTS_Kθ(hydroParam, Id_Select, K_Table, KunsatModel_Lab, N_SoilSelect::Int64, Path::String)
				println("    ~  $(Path) ~")

				N_K = Int64(length(K_Table))

			# Writting the Header
				FieldName_String =["Id", "H[mm]" ,"Kunsat[mm_s]"]
							
			# Computing K at required Ψ
				N = N_K *N_SoilSelect
            Id₂     = Vector{Int64}(undef, N)
            Ψ₂      = Vector{Float64}(undef,  N)
            Kunsat₂ = Vector{Float64}(undef, N)
            iCount  = 1
				hydroParam₂ = deepcopy(hydroParam)
				for iZ=1:N_SoilSelect
					for iK =1:N_K
						hydroParam₂.Ks[iZ] = KunsatModel_Lab[iZ]

						Id₂[iCount] = Id_Select[iZ]
						Ψ₂[iCount] = K_Table[iK]
						Kunsat₂[iCount] = kunsat.Ψ_2_KUNSAT(Ψ₂[iCount], iZ, hydroParam₂)

						iCount+=1
					end # iΨ
				end # iZ

				Output = Tables.table([string.(Id₂[1:N]) Ψ₂[1:N] Kunsat₂[1:N]])
				CSV.write(Path, Output, header=FieldName_String, delim=',')
				
		return nothing
		end  # function:  θΨ
			
	end  # module hydro
	# ............................................................


	# =============================================================
	#		MODULE: psd
	# =============================================================
	module psd
		import ...path, ...tool, ...wrc, ...param, ...cst
		import DelimitedFiles
		export PSD, θΨK_PSD

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : PSD
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function PSD(Id_Select, N_SoilSelect, paramPsd)
				println("    ~  $(path.Table_Psd) ~")

				Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(N_SoilSelect,  paramPsd)
				
				pushfirst!(FieldName_String, string("Id")) # Write the "Id" at the very begenning

				open(path.Table_Psd, "w") do io
					DelimitedFiles.write(io, [0xef,0xbb,0xbf])  # To reading utf-8 encoding in excel
					DelimitedFiles.writedlm(io,[FieldName_String] , ",",) # Header
					DelimitedFiles.writedlm(io, [Int64.(Id_Select) round.(Matrix,digits=5)], ",")
				end
				return
			end


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : θΨK_PSD
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function θΨK_PSD(hydroPsd, Id_Select, KunsatModel_Psd, N_SoilSelect)
				println("    ~  $(path.Table_θΨ_Psd) ~")

				Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(N_SoilSelect, hydroPsd)

				Matrix = hcat(Matrix, KunsatModel_Psd)

				pushfirst!(FieldName_String, string("Id")) # Write the "Id" at the very begenning
				push!(FieldName_String, "Kunsat_Model")

				Matrix =  round.(Matrix, digits=5)
				open(path.Table_θΨ_Psd, "w") do io
					DelimitedFiles.write(io, [0xef,0xbb,0xbf])  # To reading utf-8 encoding in excel
					DelimitedFiles.writedlm(io,[FieldName_String] , ",",) # Header
					DelimitedFiles.writedlm(io, [Id_Select Matrix], ",")
				end
				return
			end  # function:  θΨK_PSD



		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : PSD_θΨ_θ
		# 		Tabular values of the PSD model
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function PSD_θΨ_θ(Id_Select, N_SoilSelect, hydroPsd)
				println("    ~  $(path.Table_Psd_θΨ_θ) ~")

				N_Ψ = Int64(length(param.psd.Ψ_Table))

				# Writting the Header
					FieldName_String = Array{String}(undef, (N_Ψ))

					for i =1:N_Ψ
						FieldName_String[i] = string(param.psd.Ψ_Table[i] * cst.Mm_2_Cm) * "cm"
					end
					pushfirst!(FieldName_String, string("Id")) # Write the "Id" at the very begenning
				
				# Computing θ at required θ
					θ = Array{Float64}(undef, (N_SoilSelect, N_Ψ))
					for iZ=1:N_SoilSelect
						for iΨ =1:N_Ψ
							Ψ = param.psd.Ψ_Table[iΨ]
							θ[iZ, iΨ] = wrc.Ψ_2_θDual(Ψ, iZ, hydroPsd)
						end # iΨ
					end # iZ

				# Concatenating the 2 matrices
					θ = hcat(Id_Select, θ)
					θ = round.(θ, digits=5)

				# Writting the table
					open(path.Table_Psd_θΨ_θ, "w") do io
						DelimitedFiles.writedlm(io, [FieldName_String] , ",",) # Header
						for i = 1:length(Id_Select)
							DelimitedFiles.write(io, [0xef,0xbb,0xbf])  # To reading utf-8 encoding in excel
							DelimitedFiles.writedlm(io, [θ[i, 1:N_Ψ+1]], ",")
						end # i
					end # Path
				return
			end  # function:  θΨK_PSD

	end  # module psd
	# ............................................................
	
	

	# =============================================================
	#		MODULE: infilt
	# =============================================================
	module infilt
		import ...path, ...tool
		import DelimitedFiles
		export HYDRO_INFILT

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : HYDRO_INFILT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function HYDRO_INFILT(hydroInfilt, Id_Select, KunsatModel_Infilt, N_SoilSelect)
				println("    ~  $(path.Table_HydroInfilt) ~")

				Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(N_SoilSelect, hydroInfilt)

				Matrix = hcat(Matrix, KunsatModel_Infilt)
				
				pushfirst!(FieldName_String, string("Id")) # Write the "Id" at the very begenning
				push!(FieldName_String, string("Kunsat_θΨ"))

				Matrix =  round.(Matrix, digits=10)
				open(path.Table_HydroInfilt, "w") do io
					DelimitedFiles.write(io, [0xef,0xbb,0xbf])  # To reading utf-8 encoding in excel
					DelimitedFiles.writedlm(io,[FieldName_String] , ",",) # Header
					DelimitedFiles.writedlm(io, [Id_Select Matrix], ",")
				end
				return
			end  # function: HYDRO_INFILT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : infilt
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function INFILT(Id_Select, N_SoilSelect, infiltOutput)
				println("    ~  $(path.Table_Infilt) ~")

				Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(N_SoilSelect, infiltOutput)
				
				pushfirst!(FieldName_String, string("Id")) # Write the "Id" at the very begenning

				Matrix =  round.(Matrix, digits=5)
				open(path.Table_Infilt, "w") do io
					DelimitedFiles.write(io, [0xef,0xbb,0xbf])  # To reading utf-8 encoding in excel
					DelimitedFiles.writedlm(io,[FieldName_String] , ",",) # Header
					DelimitedFiles.writedlm(io, [Id_Select Matrix], ",")
				end
				return
			end  # function: HYDRO_INFILT
		
	end  # module: infilt
	# ...........................................................

end  # module table
# ............................................................