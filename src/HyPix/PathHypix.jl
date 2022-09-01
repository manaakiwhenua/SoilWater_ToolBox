# =============================================================
#		MODULE: pathHypix
# =============================================================
module pathHypix
	import ..option, ..sitename
	import DelimitedFiles

	export PATHHYPIX

	mutable struct PATHYPIXS
		Climate
      Dates
      Discretization
      HyPix_HydroParam
      HyPix_Param
      HyPix_VegParam
      Hydraulic_Kg
      Input_OfStep
		JulesMetadata
		ProjectName_Hypix
		SiteName_Hypix
      obsθ 

		LookUpTable_CropCoeficient
      LookUpTable_Lai

		Table_DailyClimate
		Table_Discretisation
		Table_Hydro
		Table_KΨ
		Table_Performance
		Table_Q
		Table_Signature
		Table_TimeSerie
		Table_TimeSerie_Daily
		Table_Veg
		Table_Ψ
		Table_θ
		Table_θaverage
		Table_θΨ

		Plot_Hypix_θΨK
		Plot_HypixTime
		Plot_RainfallInterception
		Plot_Se_Time
		Plot_Se_Z
		Plot_Sorptivity
		Vegetation

		Plot_OfStep
		Plot_Se_Ψ_Constrained
		Plot_θΨ_Δθ
		Plot_σ2θr
		Plot_Ψmin_Ψmax
		Plots_θ∂θ∂Ψ
	end

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : PATHHYPIX
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function PATHHYPIX(iSim)	
			# INPUT NAME OF FILE
				ProjectName_Hypix = "JULES" # "JULES"; "LYSIMETERS" 
				
				# SiteName_Hypix = "Lincoln" # "TAUPO"; "OTOROHANGA"; "WAIHOU"; "WAITOA"; "HAMILTON"; "Lincoln";
				SiteName_Hypix = sitename.SITENEAME()[iSim]
		
			# HYPIX INPUT JULES	
				JulesMetadata = "JULES_LinkingData.csv"

			# HYPIX INPUT DATA
            Climate          = "Climate_2.csv"
            Dates            = "Dates.csv"
            Discretization   = "Discretization_2.csv"
            HyPix_HydroParam = "HypixHydro.csv"
            HyPix_Param      = "HyPix_Param_2.csv"
            HyPix_VegParam   = "Vegetation.csv"
            Hydraulic_Kg     = "Hydraulic_Uni_Kg2.csv"
            Input_OfStep     = "Wof_Steps.csv"
            obsθ             = "Soilmoisture.csv"
				
			# HYPIX LOOKUPTABLE
				LookUpTable_CropCoeficient = "LookUpTable_CropCoeficient.csv"
				LookUpTable_Lai            = "LookUpTable_Lai.csv"

			# HYPIX OUTPUT TABLE
            Table_DailyClimate    = "Table_DailyClimate"
            Table_Discretisation  = "Table_Discretisation.csv"
            Table_Hydro           = "Table_Hydro"
            Table_KΨ              = "Table_KΨ"
            Table_Performance     = "Table_Performance"
            Table_Q               = "Table_Q"
            Table_Signature       = "Table_Signature"
            Table_TimeSerie       = "Table_TimeSerie"
            Table_TimeSerie_Daily = "Table_TimeSerie_Daily"
            Table_Veg             = "Table_Veg"
            Table_Ψ               = "Table_H"
            Table_θ               = "Table_Sm"
            Table_θaverage        = "Table_THETAaverage"
            Table_θΨ              = "Table_θΨ"

			# HYPIX PLOTS 
            Plot_HypixTime            = "Plot_HypixTime"
            Plot_Hypix_θΨK            = "Plot_ThetaPsiK"
            Plot_RainfallInterception = "Plot_RainfallInterception"
            Plot_Se_Time              = "Plot_Se_Time.png"
            Plot_Se_Z                 = "Plot_Se_Z.png"
            Plot_Sorptivity           = "Plot_Sorptivity"
            Vegetation                = "Plot_Vegetation"

			# HYPIX PLOT OTHERS: RESULTS
				# Plot_OfStep
            Plot_Se_Ψ_Constrained = "Plot_Se_Ψ_Constrained.svg"
            Plot_Ψmin_Ψmax        = "Plot_ΨminΨmax.svg"
            Plot_θΨ_Δθ            = "Plot_θΨ_Δθ.svg"
            Plot_σ2θr             = "Plot_θr2σ.svg"
            Plots_θ∂θ∂Ψ           = "Plot_θ∂θ∂Ψ.svg"

			# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
			# 						PROCESSING DATA
			#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				Home2 = @__DIR__

				# Moving up a level
					Home = replace(Home2,"Hypix" => "" )


				# HYPIX INPUT JULES ===
            	FileHypix_Input = Home * "//INPUT//DataHyPix//" * ProjectName_Hypix * "//"
               JulesMetadata   = FileHypix_Input * JulesMetadata

				# HYPIX INPUT DATA ===
				FileHypix_Input  = Home * "//INPUT//DataHyPix//" * ProjectName_Hypix * "//" * SiteName_Hypix * "//" * SiteName_Hypix * "_"

               Climate          = FileHypix_Input * option.hyPix.ClimateDataTimestep * "_" * Climate
               Dates            = FileHypix_Input * Dates
               Discretization   = FileHypix_Input * Discretization
               HyPix_HydroParam = FileHypix_Input * HyPix_HydroParam
               HyPix_VegParam   = FileHypix_Input * HyPix_VegParam
               Hypix_Param      = FileHypix_Input * HyPix_Param
               obsθ             = FileHypix_Input * obsθ

               Input_OfStep     = Home * "//INPUT//DataHyPix//RESULTS//"

				# HYPIX LOOKUPTABLE ===
				FileHypix_LookUpTable = Home * "//INPUT//DataHyPix//LookUpTable//"
						
               LookUpTable_CropCoeficient = FileHypix_LookUpTable * LookUpTable_CropCoeficient
               LookUpTable_Lai            = FileHypix_LookUpTable * LookUpTable_Lai

				# HYPIX OUTPUT TABLE
				FileSoilHydro_Table = Home * "//OUTPUT//Hypix//" * ProjectName_Hypix * "//" * SiteName_Hypix *"//Table//" 				
					mkpath(FileSoilHydro_Table) #Make Folder if not exist

					FileSoilHydro_Table = FileSoilHydro_Table * ProjectName_Hypix * "_"

               Table_DailyClimate    = FileSoilHydro_Table  *  SiteName_Hypix * "_"* Table_DailyClimate
               Table_Discretisation  = FileSoilHydro_Table  *  SiteName_Hypix * "_" *Table_Discretisation
               Table_Hydro           = FileSoilHydro_Table  *  SiteName_Hypix * "_" *Table_Hydro
               Table_KΨ              = FileSoilHydro_Table  *  SiteName_Hypix * "_"* Table_KΨ
               Table_Q               = FileSoilHydro_Table  *  SiteName_Hypix * "_"* Table_Q
               Table_Signature       = FileSoilHydro_Table  *  SiteName_Hypix * "_"* Table_Signature
               Table_TimeSerie       = FileSoilHydro_Table  *  SiteName_Hypix * "_"* Table_TimeSerie
               Table_TimeSerie_Daily = FileSoilHydro_Table  *  SiteName_Hypix * "_"* Table_TimeSerie_Daily
               Table_Veg             = FileSoilHydro_Table  *  SiteName_Hypix * "_"* Table_Veg
               Table_Ψ               = FileSoilHydro_Table  *  SiteName_Hypix * "_"* Table_Ψ
               Table_θ               = FileSoilHydro_Table  *  SiteName_Hypix * "_"* Table_θ
               Table_θΨ              = FileSoilHydro_Table  *  SiteName_Hypix * "_"* Table_θΨ
					
				FileSoilHydro_Table_θaverage = Home * "//OUTPUT//Hypix//" * ProjectName_Hypix * "//SoilMoistureSim//" 				
					mkpath(FileSoilHydro_Table_θaverage) #Make Folder if not exist
               Table_θaverage        = FileSoilHydro_Table_θaverage *  SiteName_Hypix * "_"* Table_θaverage

				FileSoilHydro_Table_Performace = Home * "//OUTPUT//Hypix//" * ProjectName_Hypix * "//"		
               Table_Performance     = FileSoilHydro_Table_Performace  *  SiteName_Hypix * "_"* string(iSim) * "_"* Table_Performance


				# HYPIX PLOT CORE
				FileHypix_Plot = Home * "//OUTPUT//Hypix//" * ProjectName_Hypix * "//Plots//" 
				# FileHypix_Plot = Home * "//OUTPUT//Hypix//" * ProjectName_Hypix  * "//Plots//"

					mkpath(FileHypix_Plot)

               Plot_HypixTime            = FileHypix_Plot * SiteName_Hypix  * "_" * Plot_HypixTime
               Plot_Hypix_θΨK            = FileHypix_Plot * SiteName_Hypix  * "_" * Plot_Hypix_θΨK
               Plot_RainfallInterception = FileHypix_Plot * SiteName_Hypix  * "_" * Plot_RainfallInterception
               Plot_Se_Time              = FileHypix_Plot * SiteName_Hypix  * "_" * Plot_Se_Time
               Plot_Se_Z                 = FileHypix_Plot * SiteName_Hypix  * "_" * Plot_Se_Z
               Plot_Sorptivity           = FileHypix_Plot * SiteName_Hypix  * "_" * Plot_Sorptivity
               Vegetation                = FileHypix_Plot * SiteName_Hypix  * "_" * Vegetation

				# HYPIX PLOT OTHERS: RESULTS
				FileHypix_Plot_Results = Home * "//OUTPUT//Hypix//RESULTS//"

					mkpath(FileHypix_Plot_Results)

					Plot_OfStep   = FileHypix_Plot_Results
					Plots_θ∂θ∂Ψ    = FileHypix_Plot_Results * Plots_θ∂θ∂Ψ
					Plot_Ψmin_Ψmax = FileHypix_Plot_Results * Plot_Ψmin_Ψmax
					Plot_σ2θr      = FileHypix_Plot_Results * Plot_σ2θr
					Plot_θΨ_Δθ     = FileHypix_Plot_Results * Plot_θΨ_Δθ
					Plot_Se_Ψ_Constrained = FileHypix_Plot_Results * Plot_Se_Ψ_Constrained

			# STRUCTURE
				pathHyPix = PATHYPIXS(
					Climate,
					Dates,
					Discretization,
					HyPix_HydroParam,
					HyPix_Param,
					HyPix_VegParam,
					Hydraulic_Kg,
					Input_OfStep,
					JulesMetadata,
					ProjectName_Hypix,
					SiteName_Hypix,
					obsθ, 

					LookUpTable_CropCoeficient,
					LookUpTable_Lai,

					Table_DailyClimate,
					Table_Discretisation,
					Table_Hydro,
					Table_KΨ,
					Table_Performance,
					Table_Q,
					Table_Signature,
					Table_TimeSerie,
					Table_TimeSerie_Daily,
					Table_Veg,
					Table_Ψ,
					Table_θ,
					Table_θaverage,
					Table_θΨ,

					Plot_Hypix_θΨK,
					Plot_HypixTime,
					Plot_RainfallInterception,
					Plot_Se_Time,
					Plot_Se_Z,
					Plot_Sorptivity,
					Vegetation,

					Plot_OfStep,
					Plot_Se_Ψ_Constrained,
					Plot_θΨ_Δθ,
					Plot_σ2θr,
					Plot_Ψmin_Ψmax,
					Plots_θ∂θ∂Ψ)

	return pathHyPix
	end  # function: PATHHYPIX
	
end  # module pathHyPix
# ............................................................
