# =============================================================
#		module: smap2hypix
# =============================================================
module smap2hypix
   import ..path, ..tool, ..cst, ..discretization, ..path, ..hydroStruct, ..reading, ..vegStruct, ..wrc
   import DelimitedFiles, Tables, CSV

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : SMAP_2_HYDRO
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   function SMAP_2_HYPIX(SoilName_2_SiteName,  SiteName_2_θini)

      # Reducing Ks
         Ks_FactorReduce = 0.1

      # DERIVING HYDRAULIC PARAMETER
         IgnoreSoil = "Rang_81a.2" #TODO remove 
         # Path_Output =  path.Home * "//INPUT//DataHyPix//JULES//"
         
         # Read data
            # Path_Input ="D:\\DATAraw\\JULESdata\\HydraulicParam\\Jules_HydroParam_Kosugi.csv"
             Path_Input ="D:\\DATAraw\\JULESdata\\HydraulicParam\\Jules_HydroParam_Kosugi_NoStones.csv"
            Data = DelimitedFiles.readdlm(Path_Input, ',')
            Header = Data[1,1:end]
            Data = Data[2:end,begin:end]
      
            Zlayer, N  = tool.readWrite.READ_HEADER_FAST(Data, Header, "Depth_mm")
            SoilName, ~  = tool.readWrite.READ_HEADER_FAST(Data, Header, "SoilName")
            ZrootDepth_Max, ~= tool.readWrite.READ_HEADER_FAST(Data, Header, "MaxRootingDepth_mm")

         # Determening iSite when soil changes
            iLayer_End = []
            iLayer_Start = [1]
            SoilName_Initial = SoilName[1]
            SoilName_Layer = [SoilName[1]]
            i = 1

            Nsoil = 1
            for iSoilName in SoilName
               # if soil changes
               if iSoilName ≠ SoilName_Initial
                  append!(iLayer_Start, i)
                  append!(iLayer_End, i-1)
                  push!(SoilName_Layer, iSoilName)

                  SoilName_Initial = SoilName[i] # New soil
                  Nsoil += 1
               elseif  i == N
                  append!(iLayer_End, i)  
               end  # if: name
               i += 1
            end # iSoilName
         
        
         for iSite = 1:Nsoil
            println(iSite , "=", SoilName_Layer[iSite], " =", SoilName_2_SiteName[SoilName_Layer[iSite]])
            iSiteName = SoilName_2_SiteName[SoilName_Layer[iSite]]

            Path_Output =  path.Home * "//INPUT//DataHyPix//JULES//" * iSiteName

               # SMAP-HYDRO PARAMETERS ====
                  hydroSmap = hydroStruct.HYDROSTRUCT(N) # Making a structure

                  # Abstracting data
                  hydroSmap, N_SoilSelect =  reading.READ_STRUCT(hydroSmap, Path_Input; iStart=iLayer_Start[iSite], iEnd=iLayer_End[iSite])

                  for iZ=1:N_SoilSelect
                     hydroSmap.Ks[iZ] = Ks_FactorReduce * hydroSmap.Ks[iZ]
                  end

                  Path = Path_Output * "//" * SoilName_2_SiteName[SoilName_Layer[iSite]] * "_HypixHydro.csv"

                  N_iLayers = iLayer_End[iSite] - iLayer_Start[iSite] + 1

                  TABLE_HYDRO_VEG(hydroSmap, N_iLayers, Path)

                # COMPUTING θᵢₙᵢ ===
                  Path_SmapHydro = Path_Output * "//" * iSiteName * "_HypixHydro.csv"
                     
                  Path_Output_θini =  Path_Output * "//" * iSiteName * "_ThetaIni.csv"

                  θᵢₙᵢ = SiteName_2_θini[iSiteName]
                     
                  θᵢₙᵢ_Layer = COMPUTE_θINI(hydroSmap, iSiteName, Path_Output_θini, Path_SmapHydro, θᵢₙᵢ)
 

             # DISCRETISATION ====
               # Abstracting layers per soil =====
                  Zlayer_Soil = Zlayer[iLayer_Start[iSite] : iLayer_End[iSite]]
                  Path = Path_Output *  "//" * SoilName_2_SiteName[SoilName_Layer[iSite]] * "_Layer.csv"
                  Layer = collect(1:1:length( Zlayer_Soil))

                  TABLE_DISCRETIZATION(Layer, Path, Zlayer_Soil, θᵢₙᵢ_Layer)
                  
               # Automatic Disscretizing of layers per soil =====
                  Layer, Z, θᵢₙᵢ_Cell = discretization.DISCRETISATION_AUTO(Nlayer=length(Zlayer_Soil), Zlayer=Zlayer_Soil, Zroot=800.0, θᵢₙᵢ=θᵢₙᵢ_Layer)

                  Path = Path_Output * "//" * SoilName_2_SiteName[SoilName_Layer[iSite]] * "_Discretization_2.csv"

                  TABLE_DISCRETIZATION(Layer, Path, Z, θᵢₙᵢ_Cell)


               # VEGETATION PARAMETERS
                  # Vegetation parameters per soil ====
                     vegSmap = vegStruct.VEGSTRUCT()

                     # Abstracting data
                     Path_Vegetaion ="D:\\DATAraw\\JULESdata\\Vegetation\\Vegetation.csv"

                     vegSmap.Zroot = min(vegSmap.Zroot, ZrootDepth_Max[1])

                     vegSmap, N_SoilSelect = reading.READ_STRUCT(vegSmap, Path_Vegetaion; iStart=1, iEnd=1)

                     Path = Path_Output * "//" * SoilName_2_SiteName[SoilName_Layer[iSite]] * "_Vegetation.csv"

                     TABLE_HYDRO_VEG(vegSmap, 1, Path)
            end
      
   return nothing
   end  # function: SMAP_2_HYDRO

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : COMPUTE_θINI
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function COMPUTE_θINI(hydroSmap, iSiteName, Path_Output_θini, Path_SmapHydro, θᵢₙᵢ)
         # READING HYDRAULIC PARAMETERS
            # Deriving the number of soil layers
            println(Path_SmapHydro)
               Data = DelimitedFiles.readdlm(Path_SmapHydro, ',')

               N_iZ = size(Data)[1] - 1

         println(iSiteName," ====" ,θᵢₙᵢ)

         # COMPUTING θini
            θ₁ = fill(0.0::Float64, N_iZ)

            θ₁[1] = max( min(θᵢₙᵢ, hydroSmap.θs[1]), hydroSmap.θr[1])

            Se = wrc.θ_2_Se(θ₁[1], 1, hydroSmap)

            # Assuming that all layers have the same Se
            for iZ=2:N_iZ
               θ₁[iZ] = wrc.Se_2_θ(Se, iZ, hydroSmap)
            end

         # Computing 1..N_iZ for output file
            iZ = collect(1:1:N_iZ)

         # Writing to file
            Header = ["iZ";"θini"; "Layer"]

            Output = Tables.table( [iZ θ₁ iZ])
      
            CSV.write(Path_Output_θini, Output, header=Header)	   

      return θ₁
      end  # function: COMPUTE_θINI

   
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : TABLE_DISCRETIZATION
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function TABLE_DISCRETIZATION(Layer, Path, Z, θᵢₙᵢ)
         Header = ["iZ";"Z";"θini"; "Layer"]

         iZ = collect(1:1:length(Z))

         Output = Tables.table( [iZ Z θᵢₙᵢ Layer])

         CSV.write(Path, Output, header=Header)	 
      return nothing
      end  # function: TABLE


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : TABLE_DISCRETIZATION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function TABLE_HYDRO_VEG(hydroSmap, N_iLayers, Path)
			println("			~ $(Path) ~")

			Id = 1:1:N_iLayers

			Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(N_iLayers, hydroSmap)
					
			pushfirst!(FieldName_String, string("iSite")) # Write the "Id" at the very begenning

			open(Path, "w") do io
				DelimitedFiles.write(io, [0xef,0xbb,0xbf])  # To reading utf-8 encoding in excel
				DelimitedFiles.writedlm(io,[FieldName_String] , ",",) # Header
				DelimitedFiles.writedlm(io, [Int64.(Id) Matrix], ",")
			end # open
      return nothing
		end  # function: TABLE_HYDRO

end # module smap2hypix