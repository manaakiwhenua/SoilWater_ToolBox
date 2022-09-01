# =============================================================
#		MODULE: path
# =============================================================
module path
	import ..option, ..sitename

	# NAME OF FILE
		SiteName_Soilhyro = "Smap20210226" #"VCSNSmap2"; "SFF"; "PAF"; K10KPA; Smap; Smap20210226; SmapSouthland2; CantyLysimSmap; VCSNSmap; "WaikLysim"; "Convert; "SmapNZAllSoilsSmap20210326"; "Smap20210226"
		ProjectName_Hypix = "JULES" # "JULES"; "LYSIMETERS" 
		# SiteName_Hypix = "Lincoln" # "TAUPO"; "OTOROHANGA"; "WAIHOU"; "WAITOA"; "HAMILTON"; "Lincoln";

		
		# SiteName_Hypix = sitename.SITENEAME(iSim)
		Model_Name ="A"
		Select = "SELECT_1" # Select data to model
		
	# INPUT PATH
		# Smap
			Smap                    = "Layer.csv"
			HydroParam_ThetaH       = "GUI_HydroParam.csv"

		# DATA SoilWater_Toolbox
         ConvertModel            = "TableHydro_Compiled_Homogeneous.csv"
         SmapLookupTableWettable = "LookupTable_Stone.csv"
         Id_Select               = "IdSelect.csv"
         Infiltration            = "Infiltration.csv"
         Infiltration_Param      = "Infiltration_Param.csv"
         # Kunsat                = "Kunsat_H_Ks.csv"
         Kunsat                  = "KunsatH.csv"
         Kunsat_Model            = "Kunsat_H_model.csv"
         Psd                     = "Psd.csv"
         PsdΦ                    = "PsdPorosity.csv"
         Ψθ                      = "ThetaH.csv"
         BulkDensity             = "BulkDensity.csv"
         # HydroParam_Psd          = "HydroParam_Psd.csv"
         HydroParam_Infilt       = "HydroParam_Infilt.csv"
         σ_ψM_Scenario           = "σ_ψM_Scenario.csv"

		# Temporary
         Temporary_1             = "Stratford_Min.csv"
         Temporary_2             = "Stratford_Max.csv"

		# Jules
			# JulesMetadata = "JULES_LinkingData.csv"

		# DATA HYPIX
         # Climate          = "Climate_2.csv"
         # Dates            = "Dates.csv"
         # Discretization   = "Discretization_2.csv"
         # HyPix_Param      = "HyPix_Param_2.csv"
         # HyPix_VegParam   = "Vegetation.csv"
         # HyPix_HydroParam = "HypixHydro.csv"
         Hydraulic_Kg     = "Hydraulic_Uni_Kg2.csv"
         Hydraulic_Vang   = "Hydraulic_Vang.csv"
         # Input_OfStep     = "Wof_Steps.csv"
         # obsθ             = "Soilmoisture.csv"
			
			
			if option.hyPix.θΨKmodel == :Kosugi
				Hydraulic = Hydraulic_Kg
			else option.hyPix.hypix.θΨKmodel == :vanGenuchten
				Hydraulic = Hydraulic_Vang
			end

			# DATA LOOKUPTABLE
            # LookUpTable_CropCoeficient = "LookUpTable_CropCoeficient.csv"
            # LookUpTable_Lai            = "LookUpTable_Lai.csv"

	# TABLE OUTPUT PATH
		# DATA SOIL HYDRO
         Table_HydroInfilt    = "Table_HydroInfilt.csv"
         Table_Infilt         = "Table_Infilt.csv"
         Table_Psd            = "Table_Psd.csv"
         Table_Psd_θΨ_θ       = "Table_PsdTheta.csv"
         Table_θΨK₀           = "Table_ThetaHK.csv"
         Table_θΨ_Psd         = "Table_PsdHydro.csv"
         Table_ExtraPoints_θΨ = "Table_ExtraPoints_θΨ.csv"
         Table_KosugiθΨ       = "Table_KosugiθΨ.csv"

		# HYPIX
         # Table_Discretisation  = "Table_Discretisation.csv"
         # Table_Hydro           = "Table_Hydro"
         # Table_KΨ              = "Table_KΨ"
         # Table_Performance     = "Table_Performance"
         # Table_Q               = "Table_Q"
         # Table_Signature       = "Table_Signature"
         # Table_TimeSerie       = "Table_TimeSerie"
         # Table_TimeSerie_Daily = "Table_TimeSerie_Daily"
         # Table_Veg             = "Table_Veg"
         # Table_Ψ               = "Table_H"
         # Table_θ               = "Table_Sm"
         # Table_θΨ              = "Table_θΨ"
         # Table_DailyClimate    = "Table_DailyClimate"
         # Table_θaverage        = "Table_θaverage"

	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>

	# PROCESSING
		Home = @__DIR__

		# INPUT
			# DATA SOIL HYDRO
			FileDataSoilhydro_Input = Home * "//INPUT//DataSoilHydraulic//" * SiteName_Soilhyro * "//" * SiteName_Soilhyro * "_"

			# ID_Select	
				Id_Select          = FileDataSoilhydro_Input * Id_Select

				# Smap
					Smap = FileDataSoilhydro_Input * Smap
					SmapLookupTableWettable = FileDataSoilhydro_Input * SmapLookupTableWettable

				# Convert
					ConvertModel = FileDataSoilhydro_Input * ConvertModel
					
				# BulkDensity
					BulkDensity       = FileDataSoilhydro_Input * BulkDensity
				
				#Lab
               Ψθ                = FileDataSoilhydro_Input * Ψθ
               Kunsat            = FileDataSoilhydro_Input * Kunsat
               Kunsat_Model      = FileDataSoilhydro_Input * Kunsat_Model
               HydroParam_ThetaH = FileDataSoilhydro_Input * HydroParam_ThetaH
               σ_ψM_Scenario     = FileDataSoilhydro_Input * σ_ψM_Scenario
               Temporary_1       = FileDataSoilhydro_Input * Temporary_1
               Temporary_2       = FileDataSoilhydro_Input * Temporary_2
					
				#Psd
               Psd            = FileDataSoilhydro_Input * Psd
               # HydroParam_Psd = FileDataSoilhydro_Input * HydroParam_Psd

				# Infilt
               Infiltration       = FileDataSoilhydro_Input * Infiltration
               Infiltration_Param = FileDataSoilhydro_Input * Infiltration_Param
               HydroParam_Infilt  = FileDataSoilhydro_Input * HydroParam_Infilt

			# HYPIX
				# # Jules
				# 	FileHypix_Input = Home * "//INPUT//DataHyPix//" * ProjectName_Hypix * "//" 
				# 	JulesMetadata  = FileHypix_Input * JulesMetadata

				# # Other
				# 	FileHypix_Input = Home * "//INPUT//DataHyPix//" * ProjectName_Hypix * "//" * SiteName_Hypix * "//" * SiteName_Hypix * "_" 
	
				# 	Climate        = FileHypix_Input * option.hyPix.ClimateDataTimestep * "_" * Climate
				# 	Dates          = FileHypix_Input * Dates
				# 	Discretization = FileHypix_Input * Discretization
				# 	Hydraulic      = FileHypix_Input * Hydraulic
				# 	Hypix_Param    = FileHypix_Input * HyPix_Param
					
				# 	obsθ         = FileHypix_Input * obsθ
				# 	HyPix_VegParam   =  FileHypix_Input * HyPix_VegParam
				# 	HyPix_HydroParam =  FileHypix_Input * HyPix_HydroParam
				
				# 	Input_OfStep = Home * "//INPUT//DataHyPix//RESULTS//"

			# HYPIX LOOKUPTABLE
				# FileHypix_LookUpTable = Home * "//INPUT//DataHyPix//LookUpTable//"
					
				# 	LookUpTable_CropCoeficient = FileHypix_LookUpTable * LookUpTable_CropCoeficient
				# 	LookUpTable_Lai            = FileHypix_LookUpTable * LookUpTable_Lai
		# TABLE
			# SOIL HYDRO
			FileSoilHydro_Table₁ = Home * "//OUTPUT//SoilHydro//" * SiteName_Soilhyro * "//Table//" 
			#Make Folder if not exist
			mkpath(FileSoilHydro_Table₁) 
			FileSoilHydro_Table₁ = FileSoilHydro_Table₁ * SiteName_Soilhyro * "_"

				#Lab
               Table_θΨK            = FileSoilHydro_Table₁ * string(option.hydro.HydroModel) *  "_" * string(option.hydro.σ_2_Ψm) *  "_" * Model_Name * "_" * Table_θΨK₀
               Table_θΨ_Psd         = FileSoilHydro_Table₁ * string(option.psd.HydroModel) *  "_" * string(option.hydro.σ_2_Ψm) *  "_" * Model_Name * "_" * Table_θΨ_Psd
               Table_ExtraPoints_θΨ = FileSoilHydro_Table₁ *   "_" * Table_ExtraPoints_θΨ
               Table_KosugiθΨ       = FileSoilHydro_Table₁ *   "_" * Table_KosugiθΨ

				#SMAP
					Table_Smap =  FileSoilHydro_Table₁ * "Smap.csv"
 
				#Infilt
               Table_HydroInfilt    = FileSoilHydro_Table₁ * string(option.infilt.Model) * "_" *  Model_Name  *  "_" * Table_HydroInfilt
               Table_Infilt         = FileSoilHydro_Table₁ * string(option.infilt.Model) *  "_" *  Model_Name  *  "_" *  Table_Infilt

				#Psd
					Table_Psd         = FileSoilHydro_Table₁ * string(option.psd.Model) *  "_" * Model_Name * "_" * Table_Psd
					Table_Psd_θΨ_θ    = FileSoilHydro_Table₁ * string(option.psd.HydroModel) *  "_" * Model_Name * "_" *  Table_Psd_θΨ_θ

			# HYPIX OUTPUTS
				# FileSoilHydro_Table = Home * "//OUTPUT//Hypix//" * ProjectName_Hypix * "//" * SiteName_Hypix *"//Table//" 
				
				# mkpath(FileSoilHydro_Table) #Make Folder if not exist
				# FileSoilHydro_Table = FileSoilHydro_Table * ProjectName_Hypix * "_"

            #    Table_Discretisation  = FileSoilHydro_Table  *  SiteName_Hypix * "_" *Table_Discretisation
            #    Table_Hydro           = FileSoilHydro_Table  *  SiteName_Hypix * "_" *Table_Hydro
            #    Table_Q               = FileSoilHydro_Table  *  SiteName_Hypix * "_"* Table_Q
            #    Table_Signature       = FileSoilHydro_Table  *  SiteName_Hypix * "_"* Table_Signature
            #    Table_TimeSerie       = FileSoilHydro_Table  *  SiteName_Hypix * "_"* Table_TimeSerie
            #    Table_TimeSerie_Daily = FileSoilHydro_Table  *  SiteName_Hypix * "_"* Table_TimeSerie_Daily
            #    Table_Veg             = FileSoilHydro_Table  *  SiteName_Hypix * "_"* Table_Veg
            #    Table_Performance     = FileSoilHydro_Table  *  SiteName_Hypix * "_"* Table_Performance
            #    Table_Ψ               = FileSoilHydro_Table  *  SiteName_Hypix * "_"* Table_Ψ
            #    Table_θ               = FileSoilHydro_Table  *  SiteName_Hypix * "_"* Table_θ
            #    Table_θΨ              = FileSoilHydro_Table  *  SiteName_Hypix * "_"* Table_θΨ
            #    Table_KΨ              = FileSoilHydro_Table  *  SiteName_Hypix * "_"* Table_KΨ
            #    Table_DailyClimate    = FileSoilHydro_Table  *  SiteName_Hypix * "_"* Table_DailyClimate
				# 	Table_θaverage        = FileSoilHydro_Table  *  SiteName_Hypix * "_"* Table_θaverage
			
		# PLOT
			# SOIL HYDRO
			FileSoilHydro_Plot = Home * "//OUTPUT//SoilHydro//" * SiteName_Soilhyro * "//Plots//"
				#Lab
               Plots_θΨK = FileSoilHydro_Plot * "//Lab//" 
					mkpath(Plots_θΨK)
					Plots_θΨK  = Plots_θΨK * SiteName_Soilhyro * "_"

					Plots_σΨm = FileSoilHydro_Plot * "//LabSigmaHm//" 
					mkpath(Plots_σΨm)
					Plots_σΨm  = Plots_σΨm * SiteName_Soilhyro * "_"


				#Psd
               Plot_Psd_θΨ     = FileSoilHydro_Plot * "//Psd//IMP_ThetaH//"
					mkpath(Plot_Psd_θΨ)				
               Plot_Psd_θΨ = Plot_Psd_θΨ * SiteName_Soilhyro * "_"
					
               Plots_IMP_model = FileSoilHydro_Plot * "//Psd//IMP//"
					mkpath(Plots_IMP_model)
               Plots_IMP_model = Plots_IMP_model * SiteName_Soilhyro * "_"
					
               Plots_Psd       = FileSoilHydro_Plot * "//Psd//"
					mkpath(Plots_Psd)
               Plots_Psd       = Plots_Psd * SiteName_Soilhyro * "_"

               Plots_Psd_θr    = FileSoilHydro_Plot * "//Psd//ThetaR//" 
					mkpath(Plots_Psd_θr)
					Plots_Psd_θr    = Plots_Psd_θr * "Plot_ThetaR.svg"
					
				#Infiltration					
               Plots_∑infilt_Opt        = FileSoilHydro_Plot * "//Infiltration//Optimize//"
					mkpath(Plots_∑infilt_Opt)
               Plots_∑infilt_Opt        = Plots_∑infilt_Opt * SiteName_Soilhyro * "_"
					
               Plots_∑infilt_SeIniRange = FileSoilHydro_Plot * "//Infiltration//SeIni//"
					mkpath(Plots_∑infilt_SeIniRange)
               Plots_∑infilt_SeIniRange = Plots_∑infilt_SeIniRange * SiteName_Soilhyro * "_"

               Plots_∑infilt_θΨ         = FileSoilHydro_Plot * "//Infiltration//ThetaH//"
					mkpath(Plots_∑infilt_θΨ)
               Plots_∑infilt_θΨ         = Plots_∑infilt_θΨ * SiteName_Soilhyro * "_"
					
               Plots_Sorptivity_Se      = FileSoilHydro_Plot * "//Infiltration//Sorptivity//"
					mkpath(Plots_Sorptivity_Se)
               Plots_Sorptivity_Se      = Plots_Sorptivity_Se * SiteName_Soilhyro * "_"
			
			# HYPIX
				# FileHypix_Plot = Home * "//OUTPUT//Hypix//" * ProjectName_Hypix * "//" * SiteName_Hypix *"//Plots//" 
				# # FileHypix_Plot = Home * "//OUTPUT//Hypix//" * ProjectName_Hypix  * "//Plots//"
				# mkpath(FileHypix_Plot)
				# 	Plot_HypixTime  = FileHypix_Plot * ProjectName_Hypix  * "_" * "Plot_HypixTime"
				# 	Plot_Hypix_θΨK = FileHypix_Plot * ProjectName_Hypix  * "_" * "Plot_ThetaPsiK"
				# 	Plot_Se_Time   = FileHypix_Plot * ProjectName_Hypix  * "_" * "Plot_Se_Time.png"
				# 	Plot_Se_Z      = FileHypix_Plot * ProjectName_Hypix  * "_" * "Plot_Se_Z.png"
				# 	Vegetation     = FileHypix_Plot * ProjectName_Hypix  * "_" * "Plot_Vegetation"
				# 	Plot_Sorptivity =  FileHypix_Plot * ProjectName_Hypix  * "_" * "Plot_Sorptivity"
				# 	Plot_RainfallInterception = FileHypix_Plot * ProjectName_Hypix  * "_" * "Plot_RainfallInterception"

				# FileHypix_Plot_Results = Home * "//OUTPUT//Hypix//RESULTS//"
				# mkpath(FileHypix_Plot_Results)
            #    Plot_OfStep   = FileHypix_Plot_Results
            #    Plots_θ∂θ∂Ψ    = FileHypix_Plot_Results * "Plot_θ∂θ∂Ψ.svg"
            #    Plot_Ψmin_Ψmax = FileHypix_Plot_Results * "Plot_ΨminΨmax.svg"
				# 	Plot_σ2θr      = FileHypix_Plot_Results * "Plot_θr2σ.svg"
				# 	Plot_θΨ_Δθ     = FileHypix_Plot_Results * "Plot_θΨ_Δθ.svg"
				# 	Plot_Se_Ψ_Constrained = FileHypix_Plot_Results * "Plot_Se_Ψ_Constrained.svg"

			
end  # module path
# ............................................................
