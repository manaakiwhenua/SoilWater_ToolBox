# =============================================================
#		MODULE: param
# =============================================================
	mutable struct PARAM
		globalparam
		hydro
		psd
		infilt
		hyPix 
	end
		mutable struct GLOBALPARAM
			N_iZ_Simulations
			N_iZ_Plot_Start
			N_iZ_Plot_End
		end
		mutable struct HYDROS
			Coeff_Φ_2_θs
			θs_θsMacMat
			ΨmacMat
			Ψ_Max
			Ψ_Table
			Ψ_TableComplete
			K_Table
			kg # it is a structure
			smap # it is a structureggs
		end
			mutable struct KG
				ΨmMac
				ΨmMac_Min
				ΨmMac_Max
				σMac
				∇_θsMacMat_Min
				∇_θsMacMat_Max
				Ψσ_Min
				Ψσ_Max
				Ψσ
				Pσ_1
				Pσ_2
			end
			mutable struct SMAPS
				Ψ_Table
			end
		mutable struct PSDS
			Psd_2_θr_α1
			Psd_2_θr_α1_Min
			Psd_2_θr_α1_Max
			Psd_2_θr_α2
			Psd_2_θr_α2_Min
			Psd_2_θr_α2_Max
			Psd_2_θr_Size
			Ψ_Table
			imp # it is a structure
			chang # it is a structure
		end
			mutable struct IMP
				Ψ_Max
				λ
				ξ_Max
				ξ1
				ξ1_Min
				ξ1_Max
				ξ2_Max
				∑Psd_2_ξ2_β1
				∑Psd_2_ξ2_β1_Min
				∑Psd_2_ξ2_β1_Max
				∑Psd_2_ξ2_β2
				∑Psd_2_ξ2_β2_Min
				∑Psd_2_ξ2_β2_Max
				∑Psd_2_ξ2_Size
				Subclay
				Subclay_Min
				Subclay_Max
			end
			mutable struct CHANG
				ξ1
				ξ1_Min
				ξ1_Max
			end
		mutable struct INFILTS
			SeIni_Output
			Npoint_Infilt
			ΔSlope_Err_SteadyState
		end
		mutable struct HYPIXS
			iSim_Start
			iSim_End
			Year_Start::Int64
			Month_Start::Int64
			Day_Start::Int64
			Hour_Start::Int64
			Minute_Start::Int64
			Second_Start::Int64
			Year_End::Int64
			Month_End::Int64
			Day_End::Int64
			Hour_End::Int64
			Minute_End::Int64
			Second_End::Int64
			ΔZrz_Max
			ΔZdeep_max
			Cosα
			ΔHpondMax
			Ψ_Bot
			ΔT_Min
			ΔT_Max
			N_Iter
			ΔT_Rerun
			Δθ_Max
			NewtonStepWeaken
			WaterBalanceResidual_Max
			obsθ # it is a structure
			ploting # it is a structure
		end
		mutable struct OBSΘ
			NmaxFuncEvals::Int64
			Year_Start::Int64
			Month_Start::Int64
			Day_Start::Int64
			Hour_Start::Int64
			Minute_Start::Int64
			Second_Start::Int64
			Year_End::Int64
			Month_End::Int64
			Day_End::Int64
			Hour_End::Int64
			Minute_End::Int64
			Second_End::Int64
			θobs_Uncert
		end
		mutable struct PLOT
			Year_Start
			Month_Start
			Day_Start
			Hour_Start
			Minute_Start
			Second_Start
			Year_End
			Month_End
			Day_End
			Hour_End
			Minute_End
			Second_End
			Cells_Plot
			θΨ_Table
		end
		mutable struct SIGNATURES
			Month_Winter_Start
			Month_Winter_End
		end

	# =============================================================
	#		GLOBAL PARAMETERS
	# =============================================================
		# Global param 
			N_iZ_Simulations = 1000000 # maximum number of soils to be simulated (good for testing)
			N_iZ_Plot_Start  = 1 # Starting of iZ to be plotted
			N_iZ_Plot_End    = 100000 # End of iZ to be plotted

		globalparam = GLOBALPARAM(N_iZ_Simulations, N_iZ_Plot_Start, N_iZ_Plot_End)
	

	# =============================================================
	#		HYDRO PARAMETERS
	# =============================================================
		# General parameters
			Coeff_Φ_2_θs = 0.965 # 0.95

		# Bimodal
			θs_θsMacMat = 0.05 # Bimodal if θs - θsMacMat > θs_θsMacMat

		# If σMac and ΨmMac are constants
			ΨmacMat = 100.0 # 50-150 [mm] determine when matrix and macro domain 
			Ψ_Max  = 160_000.0 # [mm] min value is 150000 mm and oven dry would be the best value 
			
		# Output in tables
			Ψ_Table = [100.0, 250000.0] # mm

			Ψ_TableComplete = [0.0, 10.0, 50.0, 1_00.0, 5_00.0, 10_00.0, 50_00.0, 100_00.0, 500_00.0, 1000_00.0, 1500_00.0,  2000_00.0] # mm

			K_Table =  [0.0, 10.0, 100.0, 500.0, 1000.0, 2000.0, 4000.0, 10000.0] # mm

		#----------------------
		# .kg KOSUGI PARAMETERS
		#----------------------
			# computing ΨmMac and σMac
				ΨmMac = √ΨmacMat  # [mm]
					ΨmMac_Min = 1.0  # 50.0[mm]
					ΨmMac_Max = ΨmacMat # 400.0[mm]

				σMac = log(ΨmacMat) / (2.0 * 3.0)

				∇_θsMacMat_Min = 0.75
				∇_θsMacMat_Max = 1.0

			# relationship between σ and Ψm
			# Parameters describing the option.σ_2_Ψm = # <:Constrained> which computes Ψm_Min, Ψm_Max from σ
				Ψσ_Min =  ΨmMac # 1.0 [mm]
				Ψσ_Max =  ΨmacMat # [mm]
				Ψσ = exp((log(Ψσ_Min) + log(Ψσ_Max)) / 2.0)  # 1.0 # [mm] used if option.σ_2_Ψm = # <:UniqueRelationship> 
				
			# Depreciated
				Pσ_1 = 1.45 # 0.5920
				Pσ_2 = 0.65 # 0.7679
	
		#----------------------
		#		smap parameters
		#----------------------
		Ψ_Table = [0.0, 500.0, 1000.0, 2000.0, 3300.0, 4000.0,  10000.0,  150000.0 ] # mm

		smap = SMAPS(Ψ_Table)

		kg = KG(ΨmMac, ΨmMac_Min, ΨmMac_Max, σMac, ∇_θsMacMat_Min, ∇_θsMacMat_Max, Ψσ_Min, Ψσ_Max, Ψσ, Pσ_1, Pσ_2)

		hydro = HYDROS(Coeff_Φ_2_θs, θs_θsMacMat, ΨmacMat, Ψ_Max, Ψ_Table, Ψ_TableComplete, K_Table, kg, smap)
	

	# =============================================================
	#		PSD module parameters
	# =============================================================
		# General parameters
			Psd_2_θr_α1    = 16.01602133125399 # α1
				Psd_2_θr_α1_Min = 0.01
				Psd_2_θr_α1_Max = 100.0
			Psd_2_θr_α2    = 2.013125380534685 # α2  
				Psd_2_θr_α2_Min = 0.001
				Psd_2_θr_α2_Max = 10.0
			Psd_2_θr_Size   = 1  # size of particle size corresponding to clay fraction
			Ψ_Table = [0.0, 500.0, 1000.0, 2000.0, 4000.0, 10000.0, 150000.0] # mm

		#----------------------
		#		IMP PARAMETERS
		#----------------------
			# Optimisation of psd
				Ψ_Max = 160000.0 # [mm] min value is 150000 mm and oven dry would be the best value  
				λ 	  = 2.0 	 # exponent of the normalised Young-Laplace capillary equation # λ = 1 for new table model 1 ######
				ξ_Max = 3.0 	 # ξ maximum physical value 

			# Intergranular mixing particle model (imp)
				ξ1 = 9.040974907360946 # 9.040974907360946 (paper)5.13729 (opt)
					ξ1_Min = 0.01 # 0.0 
					ξ1_Max =  20.0

				ξ2_Max = 0.2  

				∑Psd_2_ξ2_β1   = 0.0874077451694647 #0.0874077451694647 (paper)  0.06983 (opt) Relationship which computes ξ2 from ∑Psd
					∑Psd_2_ξ2_β1_Min = 0.001 # for new table model 4   # ξ2_Min 
					∑Psd_2_ξ2_β1_Max = 0.1   # for new table model 4   # 1.0 

				∑Psd_2_ξ2_β2   = 0.9463302042937239 #0.9463302042937239 (paper) 1.40926 (opt)
					∑Psd_2_ξ2_β2_Min = 0.1 
					∑Psd_2_ξ2_β2_Max = 5.0
					
				∑Psd_2_ξ2_Size = 2  # TODO cumulative particle size fraction corresponding to very fine silt
				Subclay        = 0.69422 #0.6934995359806453 (paper) weighted of deriving the smallest particle size
					Subclay_Min = 0.1
					Subclay_Max = 1.0
		
			#----------------------
			#		CHANG PARAMETERS
			#----------------------
				# Chang param
					ξ1 = 0.5
						ξ1_Min =  0.0 # TODO ξ1_Min chang
						ξ1_Max =  1.0 # TODO ξ1_Max chang
		
		chang = CHANG(ξ1, ξ1_Min, ξ1_Max)

		imp = IMP(Ψ_Max, λ, ξ_Max, ξ1, ξ1_Min, ξ1_Max, ξ2_Max, ∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β1_Min, ∑Psd_2_ξ2_β1_Max, ∑Psd_2_ξ2_β2, ∑Psd_2_ξ2_β2_Min, ∑Psd_2_ξ2_β2_Max, ∑Psd_2_ξ2_Size, Subclay, Subclay_Min, Subclay_Max)

		psd = PSDS(Psd_2_θr_α1, Psd_2_θr_α1_Min, Psd_2_θr_α1_Max, Psd_2_θr_α2, Psd_2_θr_α2_Min, Psd_2_θr_α2_Max, Psd_2_θr_Size, Ψ_Table, imp, chang)


	# =============================================================
	#		INFILT PARAMETERS
	# =============================================================
		# General parameters
			SeIni_Output = [0.0 0.2 0.4 0.6 0.8] # [-] Different initial Se_Ini for plotting the infiltration curves 
			Npoint_Infilt = 300 # Number of points for generating infiltration plots
			ΔSlope_Err_SteadyState = 0.5 # 0.5 Maximum error of not meeting the slope
	
		infilt = INFILTS(SeIni_Output,	Npoint_Infilt,	ΔSlope_Err_SteadyState)
	

	# =============================================================
	#		HYPIX PARAMETERS
	# =============================================================
		# Simulation wanted to run 
		  iSim_Start = 29
		  iSim_End   = 30

		# Simulations start end
         Year_Start   = 2015 # 2017
         Month_Start  = 10 #6
         Day_Start    = 18
         Hour_Start   = 9
         Minute_Start = 0
         Second_Start = 0
		   Year_End     = 2018
         Month_End    = 5
         Day_End      = 18
		   Hour_End     = 9
         Minute_End   = 0
         Second_End   = 0

		# Auto discretisation
			# If auto discretisation is selected to derive the cells sizes
			ΔZrz_Max   = 25.00 # [mm] {10-50} maximum discretisation size in the root zone
			ΔZdeep_max = 75.00 # [mm] {50-100} maximum discretisation size below the root zone

		# Slope[radiant]
		  Cosα = 1.0
		
		# Ponding
			ΔHpondMax = 15.0 # [mm] for optimisation maximum ponding depth at the end of the simulation before penalty implies 

		# Richards equation
			Ψ_Bot = 0.0 # If BottomBoundary =  [Pressure]

		# Time management
         ΔT_Min   = 10.0  # 30.0 [seconds]
         ΔT_Max   = 60.0 * 60.0 # 60.0 * 60.0 / 3.0 [seconds]

		# Iterations
         N_Iter   = 30 # 80 Maximum number of iteration before changing time step
         ΔT_Rerun = 1.2 # Allowable change of ΔT without rerun of the model > 1.0
         Δθ_Max   = 10.0^-3.3			# 3.3 1.0E-3 1.0E-4 Best result between [10^-4.5 ; 10^-5.0] smaller the smaller the time step
	
		# Newton step
			NewtonStepWeaken = 0.5 # [1 0[ when 1 then no averaging when 0 then uses Ψ from previous time k-1

		# Water balance
			WaterBalanceResidual_Max = 10.0E-11 # 10.0E-11

		#----------------------
		#	obsθ parameters
		#----------------------
			NmaxFuncEvals = 130
		
			# Calibration data start end (if available)
				Year_Start   = Year_Start
				Month_Start  = Month_Start + 1
				Day_Start    = Day_Start
				Hour_Start   = Hour_Start
				Minute_Start = Minute_Start
				Second_Start = Second_Start
				Year_End     = Year_End # Does not matter if measurements is longer than obs (corect automatically)
				Month_End    = Month_End
				Day_End      = Day_End
				Hour_End     = Hour_End
				Minute_End   = Minute_End
				Second_End   = Second_End
				θobs_Uncert = 0.03 # [cm³ cm-³]

				obsθ = OBSΘ(NmaxFuncEvals, Year_Start, Month_Start, Day_Start, Hour_Start, Minute_Start, Second_Start, Year_End, Month_End, Day_End, Hour_End, Minute_End, Second_End, θobs_Uncert)

		#----------------------
		#	 Plot \ table parameters
		#----------------------
		  Year_Start   = Year_Start
        Month_Start  = Month_Start
        Day_Start    = Day_Start
        Hour_Start   = Hour_Start
        Minute_Start = Minute_Start
        Second_Start = Second_Start
        Year_End     = Year_End
        Month_End    = Month_End
        Day_End      = Day_End
        Hour_End     = Hour_End
        Minute_End   = Minute_End
        Second_End   = Second_End
		  Cells_Plot = [1 3 5 7 9 11 13]
		  θΨ_Table = [0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 2000.0, 3000.0, 4000.0, 5000.0, 6000.0, 7000, 8000, 9000, 10000.0, 2E4, 3E4, 4E4, 5E4, 6E4, 7E4, 8E4, 9E4, 1E5, 150000.0, 2E5, 3E5, 1E6] # mm

		ploting = PLOT(Year_Start, Month_Start, Day_Start, Hour_Start, Minute_Start, Second_Start, Year_End, Month_End, Day_End, Hour_End, Minute_End, Second_End, Cells_Plot, θΨ_Table)

		#----------------------
		#		signature parameters
		#----------------------
			# Seasonality
				Month_Winter_Start = 4
				Month_Winter_End   = 10
			
	# signatures = SIGNATURES(Month_Winter_Start, Month_Winter_End)

	hyPix = HYPIXS(iSim_Start, iSim_End, Year_Start, Month_Start, Day_Start, Hour_Start, Minute_Start, Second_Start, Year_End, Month_End, Day_End, Hour_End, Minute_End, Second_End,	ΔZrz_Max, ΔZdeep_max, Cosα, ΔHpondMax, Ψ_Bot, ΔT_Min, ΔT_Max, N_Iter, ΔT_Rerun, Δθ_Max, NewtonStepWeaken, WaterBalanceResidual_Max, obsθ, ploting)


	# =============================================================
	#		PARAMETERS
	# =============================================================
		param = PARAM(globalparam, hydro, psd, infilt, hyPix)

# end param