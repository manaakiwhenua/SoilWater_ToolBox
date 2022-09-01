# =============================================================
#		module: readSmap
# =============================================================
module readSmap
   import ..path, ..tool, ..cst
   using Polynomials
   using DelimitedFiles
   export DATA2D, SMAP, ROCKFRAGMENT_WETTABLE_STRUCT

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION :SMAP
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      struct SMAP_STRUCT
         Depth        ::Vector{Float64}
         IsTopsoil    ::Vector{Int64}
         Soilname     ::Vector{String}
         RockFragment ::Vector{Float64}
         RockClass    ::Vector{String}
         RockDepth    ::Vector{Float64}
         MaxRootingDepth ::Vector{Float64}
      end
      function SMAP(Id_Select_True, N_SoilSelect)
         println("    ~  $(path.Smap) ~")

         # Read data
            Data = DelimitedFiles.readdlm(path.Smap, ',')
         # Read header
            Header = Data[1,1:end]
         # Remove first READ_ROW_SELECT
            Data = Data[2:end,begin:end]
         # Sort data
            Data = sortslices(Data, dims=1)
         
         IsTopsoil, ~  = tool.readWrite.READ_HEADER_FAST(Data, Header, "IsTopsoil")
         IsTopsoil = 	Int64.(IsTopsoil[Id_Select_True])

         Soilname, ~  = tool.readWrite.READ_HEADER_FAST(Data, Header, "Soilname")
         Soilname = Soilname[Id_Select_True] # Selecting the data

         RockClass, ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, "RockClass")
         RockClass = RockClass[Id_Select_True] # Selecting the data

         Depth, ~  = tool.readWrite.READ_HEADER_FAST(Data, Header, "depth_mm")
         Depth = Float64.(Depth[Id_Select_True]) # Selecting the data

         RockFragment, ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, "Stone_Prop")
         RockFragment = Float64.(RockFragment[Id_Select_True]) # Selecting the data

         RockDepth, ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, "RockDepth_mm")
         RockDepth = Float64.(RockDepth[Id_Select_True]) # Selecting the data

         MaxRootingDepth, ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, "MaxRootingDepth_mm")
         MaxRootingDepth = Float64.(MaxRootingDepth[Id_Select_True]) # Selecting the data

         smap = SMAP_STRUCT(Depth, IsTopsoil, Soilname, RockFragment, RockClass, RockDepth, MaxRootingDepth)			
      return smap
      end  # function: SMAP

      
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : ROCKFRAGMENT_WETTABLE
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      struct ROCKFRAGMENT_WETTABLE_STRUCT
         # RockClass::Array{String}
         RockClass_Dict::Dict{String, Int64} 
         Ψ_Rf::Array{Float64} 
         θ_Rf::Array{Float64}
         N_Ψ::Array{Int64}
         N_RockClass::Int64
         RockClass_Polynomial_Array::Array{} 
      end
      function ROCKFRAGMENT_WETTABLE()
         Path = path.SmapLookupTableWettable
         println("    ~  $(Path) ~")
         
         # Read data
            Data = DelimitedFiles.readdlm(Path, ',')
         # Read header
            Header = Data[1,1:end]
         # Remove first READ_ROW_SELECT
            Data = Data[2:end,begin:end]
         # Sort data
            # Data = sortslices(Data, dims=1, by=x->(x[1],x[2]), rev=false)

            RockClass, N_RockClass = tool.readWrite.READ_HEADER_FAST(Data, Header, "RockClass")

            RockClass_Unique = unique(RockClass)
            
            N_RockClass = length(RockClass_Unique)

            # Dictionary
            RockClass_Dict = Dict("a"=>9999)
            for i=1:N_RockClass
               RockClass_Dict[RockClass_Unique[i]] = i
            end

            Ψ₂, N₂ = tool.readWrite.READ_HEADER_FAST(Data, Header, "H[mm]")
            θ₂, ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, "Theta") 

            Ψ_Rf = zeros(Int64,(N_RockClass, 100))
            θ_Rf = zeros(Float64,(N_RockClass, 100))
            N_Ψ = zeros(Int64,(N_RockClass))

            iRockClass=1 ; iΨ=1
            for i=1:N₂
               if RockClass[i] == RockClass_Unique[iRockClass]
                  Ψ_Rf[iRockClass,iΨ] = Ψ₂[i]
                  θ_Rf[iRockClass,iΨ] = θ₂[i]
               else
                  N_Ψ[iRockClass]  = iΨ -1
                  iRockClass += 1
                  iΨ = 1
                  Ψ_Rf[iRockClass,iΨ] = Ψ₂[i]
                  θ_Rf[iRockClass,iΨ] = θ₂[i]
               end

               iΨ += 1
            end # for i=1:N

            N_Ψ[iRockClass]  = iΨ - 1

         RockClass_Polynomial_Array = []
         for iRockClass=1:N_RockClass
            RockClass_Polynomial = Polynomials.fit(log1p.(Ψ_Rf[iRockClass,1:N_Ψ[iRockClass]]), θ_Rf[iRockClass,1:N_Ψ[iRockClass]])
            X = log1p.(Ψ_Rf[iRockClass,1:N_Ψ[iRockClass]])

            Coeffs = coeffs(RockClass_Polynomial)
         
            RockClass_Polynomial_Array = push!(RockClass_Polynomial_Array, [Coeffs])
         end

      return rfWetable = ROCKFRAGMENT_WETTABLE_STRUCT(RockClass_Dict, Ψ_Rf, θ_Rf, N_Ψ, N_RockClass, RockClass_Polynomial_Array)	
      end  # function: ROCKFRAGMENT_WETTABLE


      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #		FUNCTION : DATA2D
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         function DATA2D(Path)
            println("    ~  $(Path) ~")

            # Read data
				   Data = DelimitedFiles.readdlm(Path, ',')
            # Read header
               Header = Data[1,1:end]
            # Remove first READ_ROW_SELECT
               Data = Data[2:end,begin:end]
            # Sort data
               Data = sortslices(Data, dims=1)

            # Read data of interest
               Id₂, N_SoilSelect = tool.readWrite.READ_HEADER_FAST(Data, Header, "Id")

               Soilname₂, ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, "Soilname")

               Ψdata = []
               θData = []
               for iHeader in Header
                  if occursin("wrc", iHeader)
                     θ₀, N_SoilSelect = tool.readWrite.READ_HEADER_FAST(Data, Header, iHeader)

                     iHeader = replace(iHeader, "wrc" => "")
                     iHeader = replace(iHeader, "kpa" => "")
                     iHeader = replace(iHeader, " " => "")
                     iHeader_Float=  parse(Float64, iHeader)

                     iHeader_Float = iHeader_Float * cst.kPa_2_Mm

                     append!(Ψdata, iHeader_Float)

                     try
                        θData = hcat(θData[1:N_SoilSelect, :], θ₀[1:N_SoilSelect])
                     catch
                        θData = θ₀[1:N_SoilSelect]
                     end
                  end # occursin("wrc", iHeader)
               end # for iHeader in Header

               θ_θΨ₂ = zeros(Float64, N_SoilSelect, length(Ψdata))
               Ψ_θΨ₂ = zeros(Float64, N_SoilSelect, length(Ψdata))
               N_θΨ₂ = zeros(Int64, N_SoilSelect)
    
               for iZ=1:N_SoilSelect
                  iΨ_Count = 1
                  for iΨ=1:length(Ψdata)
                     if !isnan(θData[iZ, iΨ])
                        Ψ_θΨ₂[iZ, iΨ_Count] = Ψdata[iΨ]
                        θ_θΨ₂[iZ, iΨ_Count] = θData[iZ, iΨ]
                        N_θΨ₂[iZ] += 1
                        iΨ_Count += 1
                     end #  !isnan(θData[iZ, iΨ])
                  end # iΨ
               end # iZ

         return Id₂, N_θΨ₂, Soilname₂, θ_θΨ₂, Ψ_θΨ₂
         end  # function: DATA2D
   
end  # module: readSmap
# ............................................................