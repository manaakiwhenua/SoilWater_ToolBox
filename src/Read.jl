 # =============================================================
#		MODULE: reading
# =============================================================
module reading
	import ..option, ..path, ..tool, ..param, ..vegStruct
	import  DelimitedFiles

	export ID, θΨ, KUNSATΨ, INFILTRATION, PSD, READ_STRUCT


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ID
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ID()
			println("    ~  $(path.Id_Select) ~")

			# Read data
				Data = DelimitedFiles.readdlm(path.Id_Select, ',')

				Header = Data[1,begin:end]

				Data = Data[2:end,begin:end]

				Data = sortslices(Data, dims=1)

				# Data, Header = DelimitedFiles.readdlm(path.Id_Select, ',', Any, header=true)

				Id, ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, "Id")
			
				Id_True, N_iZ_All = tool.readWrite.READ_HEADER_FAST(Data, Header, path.Select)

				Id = Int64.(Id)
			
				Id_True = Int64.(Id_True)

				Id_Select_True = convert.(Bool, Id_True)

			# Checking for errors
				for iZ=2:N_iZ_All
					if (Id[iZ] - Id[iZ-1]) < 1
						error("Id does not increase monotically at iD $(Id[iZ]) ")
					end
				end # for iZ=2:N_iZ_All
		
			N_SoilSelect = sum(Id_True)

			Id_Select = fill(0::Int64, N_SoilSelect)
			# Id_Select_True =  Array{Bool}(undef, N_iZ_All)
			iTrue = 1

			for iZ = 1:N_iZ_All
				if Id_True[iZ] == 1
					Id_Select[iTrue] = Id[iZ]
					iTrue += 1
				end	# Id_Soil == 1	
			end  # for: Id_Soil = Id_True
			
		return Id_Select, Id_Select_True, N_SoilSelect
		end  # function: ID



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : BULKDENSITY_FINEEARTH_ROCKFRAGMENT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# function BULKDENSITY_FINEEARTH_ROCKFRAGMENT(Id_Select_True, N_SoilSelect, Path)
		# 	# Load data
		# 		Data = JuliaDB.loadtable(Path)

		# 		# Getting the data
		# 			ρb_Fe_Rf = JuliaDB.select(Data, :BulkDensity_FineEarth_RockFragment_g_cm3)
		# 			ρb_Fe_Rf = ρb_Fe_Rf[Id_Select_True]

		# 			ρb_Fe = JuliaDB.select(Data,:BulkDensity_FineEarth_g_cm3)
		# 			ρb_Fe_Rf = ρb_Fe_Rf[Id_Select_True]

		# 			ρp_Fe = JuliaDB.select(Data,:ParticleDensity_FineEarth_g_cm3)
		# 			ρb_Fe_Rf = ρb_Fe_Rf[Id_Select_True]

		# 			ρ_Rf = JuliaDB.select(Data,:Density_RockFragment_g_cm3)
		# 			Rf_01 = JuliaDB.select(Data,:RockFragment_0_1)
		# 			Rf_Rlarge = JuliaDB.select(Data,:RockFragment_Rlarge_mm)
		# 			Rf_Rsmall = JuliaDB.select(Data,:RockFragment_Rsmall_mm)
		# 			Rf_Wetability = JuliaDB.select(Data, :RockFragment_Wettability_mm_mm2)

		# 	return
		# end  # function: name



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ρ_Ψθ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function BULKDENSITY(Id_Select, N_SoilSelect)
			println("    ~  $(path.BulkDensity) ~")

			# Read data
				Data = DelimitedFiles.readdlm(path.BulkDensity, ',')
			# Read header
				Header = Data[1,1:end]
			# Remove first READ_ROW_SELECT
				Data = Data[2:end,begin:end]
			# Sort data
				Data = sortslices(Data, dims=1)

			ρbSoil, ~  = tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header, "BulkDensitySoil[g_cm-3]",  N_SoilSelect)
			
			ρp_Fine, ~ = tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header, "ParticleDensity_Fine[g_cm-3]",  N_SoilSelect)

			ρ_Rock, ~  = tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header, "Density_Rock[g_cm-3]", N_SoilSelect)
			
			RockW, ~   = tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header,"Rock%", N_SoilSelect)

		return RockW, ρ_Rock, ρbSoil, ρp_Fine
		end # function: BulkDensity

	
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : θψ_FILE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function θψ_FILE(N_SoilSelect, θ_θΨ, Ψ_θΨ, N_θΨ)

			Path = path.Table_ExtraPoints_θΨ

			# Read data
				Data = DelimitedFiles.readdlm(Path, ',')
			# Read header
				Header = Data[1,1:end]
			# Remove first READ_ROW_SELECT
				Data = Data[2:end,begin:end]

				N_Ψ = Int64(length(param.hydro.Ψ_Table))

			# Writting the Header
				FieldName_String = Array{String}(undef, (N_Ψ))
				for iΨ =1:N_Ψ
					FieldName_String[iΨ] = string(Int64(param.hydro.Ψ_Table[iΨ]) ) * "mm"

					θobs, ~  = tool.readWrite.READ_HEADER_FAST(Data, Header, FieldName_String[iΨ])

					θ_θΨ =  [θobs[1:N_SoilSelect] θ_θΨ[1:N_SoilSelect,:] ]

					Ψ_Table = fill(Float64(param.hydro.Ψ_Table[iΨ]), N_SoilSelect)
				
					Ψ_θΨ = [Ψ_Table[1:N_SoilSelect] Ψ_θΨ[1:N_SoilSelect,:] ]
				end #for iΨ =1:N_Ψ

				for iZ=1:N_SoilSelect
					N_θΨ[iZ] += 2
				end
					
		return N_θΨ, θ_θΨ, Ψ_θΨ
		end  # function: θψ_FILE
	


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : θΨ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function θΨ(Id_Select, N_SoilSelect)
			println("    ~  $(path.Ψθ) ~")

			# Read data
				Data = DelimitedFiles.readdlm(path.Ψθ, ',')
			# Read header
				Header = Data[1,1:end]
			# Remove first READ_ROW_SELECT
				Data = Data[2:end,begin:end]
			# Sort data
				Data = sortslices(Data, dims=1, by=x->(x[1],x[2]), rev=false)

			# Get the data of interest
				Ψ_θΨ, N_θΨ  = tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header, "H[mm]", N_SoilSelect)
			
				θ_θΨ, ~ = tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header, "Theta", N_SoilSelect)

		return θ_θΨ, Ψ_θΨ, N_θΨ
		end  # function: θΨ


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : KUNSATΨ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function KUNSATΨ(Id_Select, N_SoilSelect)
			# Determeining where to read the data
			if isfile(path.Kunsat)
				Path = path.Kunsat
			elseif isfile(path.Kunsat_Model)
				Path = path.Kunsat_Model
			else
				error("\n SoilWater-ToolBox input error: No Kunsat data. You coud derive K(θ) from Kosugi model with option.UsePointKosugiBimodal \n")
			end
			
			println("    ~  $(Path) ~")

			# Read data
				Data = DelimitedFiles.readdlm(Path, ',')

			# Read header
				Header = Data[1,1:end]
			# Remove first READ_ROW_SELECT
				Data = Data[2:end,begin:end]
			# Sort data
				Data = sortslices(Data, dims=1, by=x->(x[1],x[2]), rev=false)
			
			Ψ_KΨ, N_KΨ = tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header, "H[mm]", N_SoilSelect)
			
         K_KΨ, ~    = tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header,"Kunsat[mm_s]", N_SoilSelect)
			
		return K_KΨ, Ψ_KΨ, N_KΨ 
		end  # function: θΨ


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : READ_STRUCT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function READ_STRUCT(structures, Path; iStart=1, iEnd=2^63 - 1)
			println("    ~  $(Path) ~")

			# Read data
				Data = DelimitedFiles.readdlm(Path, ',')
			# Read header
				Header = Data[1,1:end]
			# Remove first READ_ROW_SELECT
				Data = Data[2:end,1:end]
			# Select data of interest
				N_SoilSelect = size(Data)[1] # Initial
				iEnd= min(N_SoilSelect, iEnd)
				Data = Data[iStart:iEnd,1:end]
				N_SoilSelect = iEnd - iStart + 1 # Final

			# Reading the Model data
			for iFieldname in propertynames(structures)

				# Putting the values of Output_Vector into structures					
				try
					Output_Vector, Ndata = tool.readWrite.READ_HEADER_FAST(Data, Header, string(iFieldname))

					# Putting the values of OutPutData into the the structures hydroData
						if length(Output_Vector) == 1 # Not a vector
							Output_Vector = Output_Vector[1]
						end

						setfield!(structures, Symbol(iFieldname), Float64.(Output_Vector))
				catch
					# @warn "SoilWater-ToolBox: cannong find $iFieldname"
					Output_Vector = fill(0.0::Float64, N_SoilSelect)
					if length(Output_Vector) == 1 # Not a vector
                  Output_Vector = Output_Vector[1]
               end
	
					setfield!(structures, Symbol(iFieldname), Float64.(Output_Vector))
				end
			end

		return structures, N_SoilSelect
		end  # function: READ_STRUCT


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTRATION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		mutable struct INFILT
			RingRadius
			θ_Ini
			γ
			β
		end # struct INFILT

		function INFILTRATION(Id_Select, N_SoilSelect)
			println("    ~  $(path.Infiltration) ~")

			# Read data
				Data = DelimitedFiles.readdlm(path.Infiltration, ',')
			# Read header
				Header = Data[1,1:end]
			# Remove first READ_ROW_SELECT
				Data = Data[2:end,begin:end]
			# Sort data
				Data = sortslices(Data, dims=1)

			# Reading select data
				Tinfilt, N_Infilt = tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header, "Tinfilt[s]", N_SoilSelect)
				
				∑Infilt_Obs , ~ = tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header, "Cumul_Infiltration[mm]", N_SoilSelect)
				
			println("    ~  $(path.Infiltration_Param) ~")

			# Read data
				Data = DelimitedFiles.readdlm(path.Infiltration_Param, ',')
			# Read header
				Header = Data[1,1:end]
			# Remove first READ_ROW_SELECT
				Data = Data[2:end,begin:end]
			# Sort data
				Data = sortslices(Data, dims=1, by=x->(x[1],x[2]), rev=false)

			RingRadius , ~  = tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header, "RingRadius[mm]", N_SoilSelect; N_Point_Max=1)

			θ_Ini , ~       = tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header,"Theta_Ini[-]", N_SoilSelect; N_Point_Max=1)

			γ , ~           = tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header, "Lambda[-]", N_SoilSelect; N_Point_Max=1)

			β , ~           = tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header, "Beta[-]", N_SoilSelect; N_Point_Max=1)

			infiltParam = INFILT(RingRadius, θ_Ini, γ, β)

		return Tinfilt, ∑Infilt_Obs, N_Infilt, infiltParam
		end  # function: INFILTRATION


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PSD
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function PSD(Id_Select, N_SoilSelect) # TODO make sure that the particles are ordered from smalest to largest
			println("    ~  $(path.Psd) ~")

			# Read data
				Data = DelimitedFiles.readdlm(path.Psd, ',')
			# Read header
				Header = Data[1,1:end]
			# Remove first READ_ROW_SELECT
				Data = Data[2:end,begin:end]
			# Sort data
				Data = sortslices(Data, dims=1)

			Diameter_Psd, N_Psd 	= tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header, "Diameter[mm]",  N_SoilSelect)

			∑Psd , ~ 				= tool.readWrite.READ_ROW_SELECT(Id_Select, Data, Header,"Cumul_Psd", N_SoilSelect)

			Rpart = @. Diameter_Psd / 2.0

		return Rpart, ∑Psd, N_Psd
		end  # function: PSD



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYDRO_PARAM
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		mutable struct OPTIM
			Param_Name :: Vector{String}
			ParamOpt_Min :: Vector{Float64}
			ParamOpt_Max :: Vector{Float64}
			Param_Min :: Vector{Float64}
			Param_Max :: Vector{Float64}
			ParamOpt :: Vector{String}
			NparamOpt :: Int64
			Flag_Opt :: Bool
			ParamOpt_LogTransform :: Vector{Bool}
		end

		function HYDRO_PARAM(hydro, N_SoilSelect, Path)
		# Read data
			Data = DelimitedFiles.readdlm(Path, ',')
		# Read header
			Header = Data[1,1:end]
		# Remove first READ_ROW_SELECT
			Data = Data[2:end,begin:end]

		# Reading the Model data
			HydroModel, Ndata   = tool.readWrite.READ_HEADER_FAST(Data, Header, "MODEL")

		# Determening which parameters correspond to the selected model
		iSelectModel = [] 
		for i=1:Ndata
			if HydroModel[i] == string(option.hydro.HydroModel)
				append!(iSelectModel, i)
			end
		end

		# Reading the names of the parameters
			Param_Name, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "ParamName")
				# Selecing data
				Param_Name = Param_Name[iSelectModel]

		# Reading minimum value of the parameters
			Param_Min, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "MIN")
				# Selecing data
				Param_Min = Param_Min[iSelectModel]

		# Reading maximum value of the parameters
			Param_Max, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "MAX")
				# Selecing data
				Param_Max= Param_Max[iSelectModel]

		# Reading parameters requires log transformation [1 or 0]
			Opt_LogTransform, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "LogTransform")
				# Selecing data
				Opt_LogTransform= Opt_LogTransform[iSelectModel]

		# Reading the values of the parameters if they are not optimized
			ParamValue, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "VALUE")
				# Selecing data
				ParamValue = ParamValue[iSelectModel]

		# Reading which parameters to be optimized [1 or 0]
			Opt, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "OPT")
			# Selecing data
			Opt = Opt[iSelectModel]
			
		# Determine if we need to optimize
			if sum(Opt) ≥ 1
				Flag_Opt = true
			else
				Flag_Opt = false
			end

		# ====================================================
		ParamOpt              = []
		ParamOpt_Max          = []
		ParamOpt_Min          = []
		ParamOpt_LogTransform = []

		i = 1
		# For every hydraulic parameter
		for inParamValue in Param_Name
			# Putting the value of the parameters in hydro. Repeating the value of the parameter for all soils data: N_SoilSelect
			ParamValue_Vector = fill(Float64(ParamValue[i]), N_SoilSelect)
			setfield!(hydro, Symbol(inParamValue), ParamValue_Vector)

			# θsMacMat value depends on θs
			if  Symbol(inParamValue) == :θsMacMat_ƞ
				for iZ = 1:N_SoilSelect 
					hydro.θsMacMat[iZ] =  hydro.θs[iZ] * hydro.θsMacMat_ƞ[iZ]
				end
			end # Symbol(inParamValue) == :θsMacMat_ƞ

			# Putting the minimum value in the parameter
				ParamValue_Vector = fill(Float64(Param_Min[i]), N_SoilSelect)
				setfield!(hydro, Symbol(inParamValue * "_Min"), ParamValue_Vector)

			# Putting the maximum value in the parameter
				ParamValue_Vector = fill(Float64(Param_Max[i]), N_SoilSelect)
				setfield!(hydro, Symbol(inParamValue * "_Max"), ParamValue_Vector)
	
			# ParamValue to optimize. The complication is that there are different layers of hydraulic parameters which can be optimized.  
			if Opt[i] == 1

				# appending the values of the parameters
				append!(ParamOpt, [Param_Name[i]])

				append!(ParamOpt_Min, Param_Min[i])
				
				append!(ParamOpt_Max, Param_Max[i])

				# Appending name of param to perform logTransform if optimized
				if Opt_LogTransform[i] == 1
					append!(ParamOpt_LogTransform, [true])
				else
					append!(ParamOpt_LogTransform, [false])
				end

				if Param_Min[i] > Param_Max[i]
					error("LabOpt ERROR: $(Param_Min[i]) < $(ParamValue[i]) < $(Param_Max[i]) !")
				end
			end # if Flag_Opt

			i += 1
		end # for loop

		# Number of parameters to be optimised
			NparamOpt = length(ParamOpt)
	
		# Putting all the in mutable structure
			optim = OPTIM(Param_Name,ParamOpt_Min,ParamOpt_Max,Param_Min,Param_Max,ParamOpt,NparamOpt,Flag_Opt,ParamOpt_LogTransform)

		if Flag_Opt == true
			println("		=== === Optimizing the following parameters === ===")
			println("			Model=" , option.hydro.HydroModel)
			println("			NparamOpt=" , NparamOpt)
			println("			ParamOpt= " ,  optim.ParamOpt)
			println("			Min_Value= " , optim.ParamOpt_Min)
			println("			Max_Value= " , optim.ParamOpt_Max)
			println("			LogTransform = " , optim.ParamOpt_LogTransform)
			println("		=== === ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ === === \n \n")
		end

	return hydro, optim
	end  # function: HydroParam_ThetaH

end  # module: reading
# ............................................................		