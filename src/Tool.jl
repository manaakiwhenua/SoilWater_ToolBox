# =============================================================
#		MODULE: tool
# =============================================================
module tool
	# =============================================================
	#		module: normalize
	# =============================================================
	module norm
	
		export ∇NORM_2_PARAMETER

			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			#		FUNCTION : HYDRO_ADJEUSTMENTS
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function ∇NORM_2_PARAMETER(∇P, P_Min, P_Max)
				return P = ∇P * (P_Max - P_Min) + P_Min
			end  # function: HYDRO_ADJEUSTMENTS
		
	end  # module: normalize
	# ............................................................

	# =============================================================
	#		MODULE: array
	# =============================================================
	module array
		export SEARCH_INDEX, SEARCH_INDEX2

		function SEARCH_INDEX(Array, SearchValue)
			N = length(Array)
			iSearchValue = 1
			Value_SearchValue=1.
			Err_2 = 100000000000000.
			
			for i = 1:N
				Err_1 = abs(Array[i] - SearchValue)
			
				if Err_1 < Err_2
					iSearchValue = i
					Value_SearchValue = Array[i]
					Err_2 = Err_1
				end
			end
			return iSearchValue
		end # function SEARCH_INDEX

		function SEARCH_INDEX2(Find, Array, N )
			i = 2
			FlagBreak = false
			while !(FlagBreak)
				if (Array[i-1] <= Find <= Array[i]) || (i == N) 
					FlagBreak = true
					break
				else 
					i += 1
					FlagBreak = false
				end # if
			end # while
			
			return i, Array[i] 
		end  # function SEARCH_INDEX2

	end  # module: array
	# ............................................................


	# =============================================================
	#		MODULE: readWrite
	# =============================================================
	module readWrite
		import ...param, ..tool
		import DelimitedFiles
		import Dates: value, DateTime
		export FIELDNAME_2_STRUCT_VECT, STRUCT_2_FIELDNAME, READ_HEADER, READ_ROW_SELECT, DATA_2_ΔTnew, STRUCT_2_FIELDNAME_PARAM

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : READ_HEADER_FAST
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function READ_HEADER_FAST(Data, Header, Name)
			N_X, N_Y = size(Data) # Size of the array

			Header = reshape(Header, N_Y, 1)
			
		# Getting the column which matches the name of the header
			Name = replace(Name, " " => "") # Remove white spaces
			try
				iColumn = Int64(findfirst(isequal(Name), Header)[1])
				global Data_Output = Data[1:N_X,iColumn]
			catch
				println(Header)
				error("\n          SOILWATERTOOLBOX ERROR: cannot find  $Name   \n \n")
			end
		return Data_Output, N_X
		end # function READ_HEADER

		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : READ_HEADER
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function READ_HEADER(Path, Name)
				# Read data
					Data =  DelimitedFiles.readdlm(Path, ',')
					N_X, N_Y = size(Data) # Size of the array
					
				# Reading header
					Header = fill("", N_Y)
					for i = 1:N_Y
						Header[i] = Data[1,i]
					end

				# Getting the column which matches the name of the header
					Name = replace(Name, " " => "") # Remove white spaces
	
					try
						global Data_Output = Data[2:N_X,findfirst(isequal(Name), Header)]
					catch
						println(Header)
						error("\n \n SOILWATERTOOLBOX ERROR: cannot find  $Name  in $Path \n \n")
					end

					N_X -= 1 # To take consideration of the header

			return Data_Output, N_X
			end # function READ_HEADER


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : READ_ROW_SELECT
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function READ_ROW_SELECT( Id_Select::Vector{Int64}, Data, Header,  Name::String, N_SoilSelect::Int64; N_Point_Max=1000, DataSelect=true)
			# READ DATA
				N_X, N_Y = size(Data) # Size of the array

			# Get the ID of the data
				Id_Data = Int64.(Data[1:end,1])

			# Getting the column which matches the name of the header
				Name = replace(Name, " " => "") # Remove white spaces

				iColumn = Int64(findfirst(isequal(Name), Header)[1])
		
				Data_Output = Float64.(Data[1:end,iColumn])

			# ===========================================
			# Only keeping data which is selected
			# ===========================================
				if DataSelect == true
					Data_Select = Array{Float64}(undef, (N_SoilSelect, N_Point_Max))
					N_Point = zeros(Int64, N_SoilSelect) 
					
					iSelect = 1; iPoint = 1
					# For all soils in the file
					for i = 1:N_X
						if Id_Data[i] == Id_Select[iSelect] # Only append Ids which correspond to the selected one
							Data_Select[iSelect,iPoint] = Data_Output[i]
							# append!(Data_Select, Data_Output[i])
							N_Point[iSelect] += 1
							iPoint += 1
						end

						# Since there are many Data_Output with the same Id only update Id_Select if we are changing soils and Id_Select[iSelect] == Id_Data[i]
						if i ≤ N_X -1
							if Id_Data[i+1] > Id_Data[i] && Id_Select[iSelect] == Id_Data[i] && iSelect ≤ N_SoilSelect -1
								iSelect += 1
								iPoint = 1
							end # if:
						end # if: i â‰¤ N_X
					end # for: i = 1:N_X
				else
					Data_Select = Data_Output
					N_Point = N_X
				end # DataSelect

		return Data_Select, N_Point
		end # function READ_ROW_SELECT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : FIELDNAME_2_STRUC
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function FIELDNAME_2_STRUCT_VECT(Structure, NameStruct)
				N_FieldName = length(fieldnames(Structure))

				FieldName_String = Array{Symbol}(undef, (N_FieldName))
				i = 1
				for FieldNames in fieldnames(Structure)
					FieldName_String[i] = FieldNames 
					i += 1
				end

				NameStruct.FieldName = FieldName_String

			return NameStruct
			end  # function: FIELDNAMES


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : STRUCT_2_FIELDNAMES
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function STRUCT_2_FIELDNAME(N_SoilSelect, Structure)

				FieldName_Array = propertynames(Structure)

				N_FieldName = length(FieldName_Array)

				# Matrix
					Matrix = Array{Float64}(undef, (N_SoilSelect, N_FieldName))

					i = 1
					for FieldName in FieldName_Array
						Struct_Array = getfield(Structure, FieldName)
		
						Matrix[1:N_SoilSelect,i] .= Struct_Array
						i += 1
					end
				
				# HEADER
					FieldName_String = Array{String}(undef, N_FieldName)
					i=1
					for FieldNames in FieldName_Array
						FieldName_String[i] =  String(FieldNames)
						i += 1
					end

				return Matrix, FieldName_String
				end # function STRUCT_2_FIELDNAME



		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : STRUCT_2_FIELDNAMES
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function STRUCT_2_FIELDNAME_PARAM(Structure)
			N_FieldName = length(Structure.FieldName) - 1

			Matrix = Array{Float64}(undef, (N_FieldName))
			
			for i=1:N_FieldName
				Struct_Array = getfield(Structure, Structure.FieldName[i])
				Matrix[i] = Struct_Array
			end

			FieldName_String = Array{String}(undef, N_FieldName)
			i=1
			for FieldNames in Structure.FieldName
				FieldName_String[i] =  String(FieldNames)
				if i == N_FieldName
					break
				end
				i += 1
			end
			return Matrix, FieldName_String
		end # function STRUCT_2_FIELDNAME
			
			
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : matching data with different timesteps
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function DATA_2_ΔTnew(∑T, N_iT, ∑T_Reduced)
			True = falses(N_iT) # Reserving memory
			iT = 1
			for iTreduced = 1:length(∑T_Reduced)	
				while ∑T_Reduced[iTreduced] ≥ ∑T[iT]
					iT += 1
				end

				# @show iT
				if iT ≤ N_iT
					True[iT] = true
				else
					exit
				end
			end # for iTθ

			# More accurate ∑T_Reduced
			∑T_Reduced = ∑T[True[1:N_iT]]

			N_∑T_Reduced = count(True[1:N_iT])

			return N_∑T_Reduced, True
		end  # function: function SELECT_OUTPUT_ΔT
		
	end  # module readWrite ************************
	# ............................................................
end  # module tool
# ............................................................