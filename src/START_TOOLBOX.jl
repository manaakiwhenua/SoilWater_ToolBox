##========================================================================================
##                                                                                      ##
##                                 Soil Water ToolBox                                   ##
##                                                                                     ##
##========================================================================================

using Suppressor
using Revise

@suppress begin
	include("Option.jl")
	# Install packages to run program
	if option.globalopt.DownloadPackage
		include("Packages.jl")
	end # option.globalopt.DownloadPackage
	include("Tool.jl")
	include("Hypix\\Other\\SiteName.jl")
	include("Path.jl")
	include("Cst.jl")
	include("Param.jl")
	include("Hydro\\ΨminΨmax.jl")
	include("Hydro\\HydroStruct.jl")
	include("Hydro\\HydroRelation.jl")
	include("Hydro\\WaterRetentionCurve.jl")
	include("Optim\\Optimize.jl")
	include("Read.jl")
	if !(option.globalopt.Hypix)
		include("Hydro\\TotalPorosity.jl")
	end
	include("Checking.jl")
	include("Hydro\\Kunsat.jl")
	include("Stats.jl")
	if !(option.globalopt.Hypix)
		include("Table.jl")
		include("Psd\\PsdThetar.jl")
	end

	if option.globalopt.Smap
		include("Smap\\StoneSmap.jl")
		include("Smap\\ReadSmap.jl")
		include("Smap\\PlotSmap.jl")
		include("Smap\\TableSmap.jl")
	end

	if option.globalopt.θΨ ≠ :No && option.globalopt.θΨ ≠ :File &&  !(option.globalopt.Hypix)
		# include("HydroLab\\HydrolabInitialize.jl")	
		include("HydroLab\\ObjectiveFunction_Lab.jl")
		# include("HydroLab\\START_Lab.jl")
		include("HydroLab\\HydrolabOpt.jl")
		include("Hypix\\Other\\PlotOther.jl")
	end
	
	if option.globalopt.Infilt
		if option.infilt.SortivityVersion == :NonInfinity
			include("Infilt\\SorptivityNonInfinity.jl")
		elseif option.infilt.SortivityVersion == :Traditional
			include("Infilt\\SorptivityTraditional.jl")
		end
		
		include("Infilt\\BestFunc.jl")
		include("Infilt\\ObjectiveFunction_Best.jl")
		include("Infilt\\QuasiExact.jl")
		include("Infilt\\TimeTransSteady.jl")
		include("Infilt\\InfiltStruct.jl")
		include("Infilt\\InfiltInitialize.jl")
		include("Infilt\\START_Infilt.jl")
	end # option.globalopt.Infilt
	
	if option.globalopt.Psd
		include("Psd\\PsdStruct.jl")
		include("Psd\\PsdInitialize.jl")
		include("Psd\\PsdFunc.jl")
		include("Psd\\PsdOpt.jl")
		include("Psd\\START_PSD.jl")
	end # option.globalopt.Psd

	if option.globalopt.Ploting && !(option.globalopt.Hypix)
		include("Plot.jl")
	end # option.globalopt.Ploting

	println(option.globalopt.Hypix)
	if option.globalopt.Hypix
		include("Hypix\\PathHypix.jl")
		include("Hydro\\ΨminΨmax.jl")
		include("Infilt\\SorptivityNonInfinity.jl")
		include("Hypix\\Interpolate.jl")
		include("Hypix\\Opt\\ThetaObs.jl")
		include("HyPix\\ThetaIni.jl")
		# include("Hypix\\Opt\\Signature.jl")
		include("Hypix\\Opt\\OfHypix.jl")
		include("Hypix\\TableHypix.jl")
		include("Hypix\\VegStruct.jl")
		include("Read.jl")
		include("Hypix\\HorizonLayer.jl")
		include("Hypix\\ReadHypix.jl")
		include("Hypix\\RainfalIntercept.jl")
		include("Hypix\\Flux.jl")
		include("Hypix\\Discretisation.jl")
		include("Hypix\\Ponding.jl")
		include("Hypix\\Residual.jl")
		include("Hypix\\ChangeOutputTimeStep.jl")
		include("Hypix\\TimeStep.jl")
		include("Hypix\\Richard.jl")
		include("Hypix\\WaterBalance.jl")
		include("Hypix\\Evaporation.jl")
		include("Hypix\\RootWaterUptake.jl")
		include("Hypix\\CheckError.jl")
		include("Hypix\\Pet.jl")
		include("HyPix\\Other\\ThetaAverage.jl")
		include("Hypix\\Memory.jl")
		include("Hypix\\Climate.jl")
		if option.globalopt.Ploting
			include("Hypix\\Other\\PlotOther.jl")
			include("Hypix\\PlotHypix.jl")
		end
		include("Hypix\\HypixModel.jl")
		include("Hypix\\Opt\\HypixOpt.jl")
		include("Hypix\\Hypix_Start.jl")
	end  # if: option.globalopt.Hypix

	if option.globalopt.Jules
		include("Hypix\\PathHypix.jl")
		include("Hypix\\VegStruct.jl")
		include("Hypix\\Discretisation.jl")
		include("Jules\\Jules.jl")
		include("HyPix\\ThetaIni.jl")
		include("Smap\\Smap2Hypix.jl")		
	end  # if: option.Temporay
end # Suppressor 


# ===============================================================
#		FUNCTION : START_TOOLBOX
# ==============================================================
function START_TOOLBOX()

	# READING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		if option.globalopt.HydroTranslateModel

			# Creating 
			hydroTranslate = hydroStruct.HYDROSTRUCT(1000)
			
			hydroTranslate, N_SoilSelect = reading.READ_STRUCT(hydroTranslate, path.ConvertModel)
			
			# Temporary Id
				Id_Select = collect(1:1:N_SoilSelect)
		
			# Deriving a table of θ(Ψ)
				table.hydroLab.TABLE_EXTRAPOINTS_θΨ(hydroTranslate, Id_Select, N_SoilSelect, path.Ψθ, param.hydro.Ψ_TableComplete; Orientation="Vertical")
			
			# Deriving a table of K(θ)
				table.hydroLab.TABLE_EXTRAPOINTS_Kθ(hydroTranslate, Id_Select, param.hydro.Ψ_TableComplete, hydroTranslate.Ks[1:N_SoilSelect], N_SoilSelect::Int64, path.Kunsat)

			# Creating an Id output required by the program
				table.TABLE_ID(N_SoilSelect::Int64, path.Id_Select)
			
		elseif !(option.globalopt.Hypix)
			# Selecting soils of interest
				Id_Select, Id_Select_True, N_SoilSelect = reading.ID()

			# Determine the soils to simulale
				N_SoilSelect = Int(min(N_SoilSelect, param.globalparam.N_iZ_Simulations))
				Id_Select = Id_Select[1:N_SoilSelect]

			# Reading bulk density
				if option.globalopt.BulkDensity
					RockW, ρ_Rock, ρbSoil, ρp_Fine = reading.BULKDENSITY(Id_Select, N_SoilSelect)
				end

			# Reading θ(Ψ)
				if option.globalopt.θΨ ≠ :No # <> = <> = <> = <> = <> = <>
					θ_θΨ, Ψ_θΨ, N_θΨ = reading.θΨ(Id_Select, N_SoilSelect)
				else
					θ_θΨ = [] 
					Ψ_θΨ = [] 
					N_θΨ = 0.0
				end
		
			# Reading K(θ)
				if option.hydro.KunsatΨ # <>=<>=<>=<>=<>
					K_KΨ, Ψ_KΨ, N_KΨ = reading.KUNSATΨ(Id_Select, N_SoilSelect)
				else
					Ψ_KΨ = []
					K_KΨ = []
					N_KΨ = 0
				end # option.hydro.KunsatΨ

			# Reading Particle Size distribution
				if option.globalopt.Psd
					Rpart, ∑Psd, N_Psd = reading.PSD(Id_Select, N_SoilSelect)
				else
					∑Psd = zeros(Float64, N_SoilSelect, 1)
					Rpart = zeros(Float64, N_SoilSelect, 1)
					N_Psd = zeros(Float64, N_SoilSelect)	
				end # option.globalopt.Psd
		
			# Reading infiltration	
				if option.globalopt.Infilt # <>=<>=<>=<>=<>
					Tinfilt, ∑Infilt_Obs, N_Infilt, infiltParam  = reading.INFILTRATION(Id_Select, N_SoilSelect)
				end # option.globalopt.Infilt

			# Reading Smap data
				if option.globalopt.Smap
					smap = readSmap.SMAP(Id_Select_True, N_SoilSelect)
					rfWetable = readSmap.ROCKFRAGMENT_WETTABLE()
				end # option.globalopt.Smap

			else # TODO: Needs to be removed
				N_SoilSelect = 1
				
			end # Option
			

		# END READING ................................................................
	if option.globalopt.Smap
		if option.smap.CorrectStone
			println("=\n  				~~~~ Stone correction ~~~~~ \n")
			θ_θΨ = stoneCorrection.STONECORRECTION(N_SoilSelect, N_θΨ, smap, θ_θΨ, Ψ_θΨ)
		end
		if option.smap.CorrectStoneWetability
			θ_θΨ = stoneCorrection.STONECORRECTION_WETTABLE(N_SoilSelect, N_θΨ, rfWetable, smap, θ_θΨ, Ψ_θΨ)
		end

		if option.smap.UsePointKosugiBimodal
			N_θΨ, θ_θΨ, Ψ_θΨ = reading.θψ_FILE(N_SoilSelect, θ_θΨ, Ψ_θΨ, N_θΨ)
		end
	end

	if option.globalopt.Jules
		SoilName_2_SiteName,  SiteName_2_θini = jules.START_JULES()
		smap2hypix.SMAP_2_HYPIX(SoilName_2_SiteName,  SiteName_2_θini)
		
	end  # if: option.START_JULES()

	if option.globalopt.θΨ ≠ :No # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	println("\n === START: DERIVING HYDRO PARAMETERS  === \n")
	println("         ===== Model Name= $(path.Model_Name) =====")
		
		# INITIALIZES HYDRAULIC PARAMETERS STRUCT INDEPENDENTLY OF THE SELECTED MODEL
			hydro = hydroStruct.HYDROSTRUCT(N_SoilSelect)

			hydroOther = hydroStruct.HYDRO_OTHERS(N_SoilSelect)

			hydro, optim = reading.HYDRO_PARAM(hydro, N_SoilSelect, path.HydroParam_ThetaH)

			checking.CHECKING(optim)

			if option.globalopt.Smap
				hydroParam, optim = stoneCorrection.STONECORRECTION_HYDRO(hydro, N_SoilSelect, optim, smap)
			end

			# plotOther.PLOT_σ_2_θr()
			# plotOther.SE_Ψ_CONSTRAINED()
			# plotOther.σ_ψM_SCEARIO()

		# If the hydraulic parameters were already derived than get the data from file instead or rerunning the model	
		if option.globalopt.θΨ == :File
			println("    ~ HydroLab HydroParam reading from file ~")
			hydro = reading.HYDROPARAM(Id_Select, N_SoilSelect, hydro)
		else
			# Total Porosity= Φ
				if option.globalopt.BulkDensity
					hydro.Φ = Φ.ρB_2_Φ(N_SoilSelect, RockW, ρ_Rock, ρbSoil, ρp_Fine)
				end

			if option.hydro.KunsatΨ
				hydro, hydroOther = hydrolabOpt.HYDROLABOPT_START(N_SoilSelect=N_SoilSelect, ∑Psd=∑Psd, θ_θΨ=θ_θΨ, Ψ_θΨ=Ψ_θΨ, N_θΨ=N_θΨ, K_KΨ=K_KΨ, Ψ_KΨ=Ψ_KΨ, N_KΨ=N_KΨ, hydro=hydro, hydroOther=hydroOther, optionHydro=option.hydro, optim=optim)

			else
				hydro, hydroOther =  hydrolabOpt.HYDROLABOPT_START(N_SoilSelect=N_SoilSelect, ∑Psd=∑Psd, θ_θΨ=θ_θΨ, Ψ_θΨ=Ψ_θΨ, N_θΨ=N_θΨ, hydro=hydro, hydroOther=hydroOther, optionHydro=option.hydro, optim=optim)
			end # option.hydro.KunsatΨ

			# SPATIAL CASE FOR BROOKS AND COREY
				if option.hydro.HydroModel==:BrooksCorey || option.hydro.HydroModel==:ClappHornberger
					for iZ=1:N_SoilSelect
						hydro.Ψga[iZ] = wrc.GREEN_AMPT(iZ, hydro)
					end
				end #  option.hydro.HydroModel
		end # option.globalopt.θΨ == :File


		# Deriving Kunsat from θ(Ψ)
		"""Pollacco, J.A.P., Webb, T., McNeill, S., Hu, W., Carrick, S., Hewitt, A., Lilburne, L., 2017. Saturated hydraulic conductivity model computed from bimodal water retention curves for a range of New Zealand soils. Hydrol. Earth Syst. Sci. 21, 2725–2737. https://doi.org/10.5194/hess-21-2725-2017"""
			KunsatModel_Lab = fill(0.0::Float64, N_SoilSelect)
			if option.hydro.HydroModel==:Kosugi 
			println("	=== START: Dering Ks from lab θ(Ψ) data ===")
				for iZ=1:N_SoilSelect
					if option.globalopt.Smap
						if smap.IsTopsoil[iZ] == 1
							TopsoilSubsoil="Topsoil"
						else
							TopsoilSubsoil="Subsoil"
						end
					else
						TopsoilSubsoil="Topsoil"
					end
					KunsatModel_Lab[iZ] = kunsat.θΨ_2_KUNSAT(0.9999, iZ, hydro, 0.0; TopsoilSubsoil=TopsoilSubsoil)

					if option.hydro.KunsatΨ == false
						hydro.Ks[iZ] = KunsatModel_Lab[iZ]
					end
				end # for
			println("		~~~ END: Dering Ks from lab θ(Ψ) data ~~~ \n")
			end

	println("=== END  : DERIVING HYDRO PARAMETERS  === \n")
	else
		hydro = []
	end


	if option.globalopt.Psd  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	println("=== START: PSD MODEL  ===")
		# Structure of hydroPsd
			hydroPsd = hydroStruct.HYDROSTRUCT(N_SoilSelect)
			hydroOther_Psd = hydroStruct.HYDRO_OTHERS(N_SoilSelect)
			hydroPsd, optim_Psd = reading.HYDRO_PARAM(hydroPsd, N_SoilSelect, path.HydroParam_ThetaH)

		# Total Porosity= Φ
		if option.globalopt.BulkDensity
			hydroPsd.Φ = Φ.ρB_2_Φ(N_SoilSelect, RockW, ρ_Rock, ρbSoil, ρp_Fine)
		end

		# PSD model
			paramPsd, N_Psd, θ_Rpart, Ψ_Rpart, Psd, hydroPsd = psd.START_PSD(∑Psd, hydro, hydroPsd, N_Psd, N_SoilSelect, N_θΨ, Rpart, θ_θΨ, Ψ_θΨ)

		KunsatModel_Psd = fill(0.0::Float64, N_SoilSelect)

		if  option.psd.HydroParam
			hydroPsd, hydroOther_Psd = hydrolabOpt.HYDROLABOPT_START(N_SoilSelect=N_SoilSelect, ∑Psd=∑Psd, θ_θΨ=θ_Rpart, Ψ_θΨ=Ψ_Rpart, N_θΨ=N_Psd, hydro=hydroPsd, hydroOther=hydroOther_Psd, optionHydro=option.psd, optim=optim_Psd)
		end

	
		println("=== END  : DERIVING HYDRO PARAMETERS  === \n")

	println("=== END : PSD MODEL  === \n")
	else
		θ_Rpart = zeros(Float64, N_SoilSelect,1)
		Ψ_Rpart = zeros(Float64, N_SoilSelect,1)
		hydroPsd = hydroStruct.HYDROSTRUCT(N_SoilSelect)
		N_Psd = zeros(Float64, N_SoilSelect)

	end # option.globalopt.Psd ...............................................................................

	
	if option.globalopt.Infilt  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	println("=== START: INFILTRATION  ===")
		# Structure of hydroInfilt
			hydroInfilt = hydroStruct.HYDROSTRUCT(N_SoilSelect)
			hydroOther_Infilt = hydroStruct.HYDRO_OTHERS(N_SoilSelect)
			hydroInfilt, optim_Infilt = reading.HYDRO_PARAM(hydroPsd, N_SoilSelect, path.HydroParam_Infilt)

		# Total Porosity= Φ
			hydroInfilt.Φ = Φ.ρB_2_Φ(N_SoilSelect, RockW, ρ_Rock, ρbSoil, ρp_Fine)

		# Running infiltration model
			infiltOutput, hydroInfilt, ∑Infilt_3D, ∑Infilt_1D = infilt.START_INFILTRATION(∑Infilt_Obs, ∑Psd, hydro, hydroInfilt, Id_Select, infiltParam, N_Infilt, N_SoilSelect, Tinfilt)

	println("=== END  : INFILTRATION  === \n")
	else
		hydroInfilt = []
	end # option.globalopt.Infilt

	if option.globalopt.Hypix
		hypix.HYPIX_START()
	end # option.globalopt.Hypix


	# TABLES OUTPUT ======================================================================================
		if option.globalopt.θΨ ≠ :No && option.globalopt.θΨ ≠ :File # <>=<>=<>=<>=<>

			if !(option.globalopt.Smap)
				table.hydroLab.θΨK(hydro, hydroOther, Id_Select[1:N_SoilSelect], KunsatModel_Lab, N_SoilSelect)
			else
				tableSmap.θΨK(hydro, hydroOther, Id_Select[1:N_SoilSelect], KunsatModel_Lab, N_SoilSelect, smap)

				if option.smap.AddPointKosugiBimodal && option.hydro.HydroModel == :Kosugi && option.hydro.σ_2_Ψm == :Constrained
					# Extra points in θ(Ψ) to reduce none uniqueness
					table.hydroLab.TABLE_EXTRAPOINTS_θΨ(hydro, Id_Select, N_SoilSelect, path.Table_ExtraPoints_θΨ, param.hydro.Ψ_Table)
	
					# Extra points required by TopNet
						table.hydroLab.TABLE_EXTRAPOINTS_θΨ(hydro, Id_Select, N_SoilSelect, path.Table_KosugiθΨ, param.hydro.smap.Ψ_Table)

						table.hydroLab.TABLE_EXTRAPOINTS_Kθ(hydro, Id_Select, param.hydro.K_Table, KunsatModel_Lab, N_SoilSelect, path.Kunsat_Model)
				end

				if option.smap.CombineData
					tableSmap.SMAP(Id_Select, N_SoilSelect, smap)
				end
			end

		end # option.globalopt.θΨ 

		if option.globalopt.Psd # <>=<>=<>=<>=<>
			table.psd.PSD(Id_Select[1:N_SoilSelect], N_SoilSelect, paramPsd)

			if option.psd.HydroParam  && option.psd.HydroParam
				table.psd.θΨK_PSD(hydroPsd, Id_Select, KunsatModel_Psd, N_SoilSelect)
			end
			
			if option.psd.Table_Psd_θΨ_θ
				table.psd.PSD_θΨ_θ(Id_Select, N_SoilSelect, hydroPsd)
			end
		end # option.globalopt.Psd

		if option.globalopt.Infilt # <>=<>=<>=<>=<>
			table.infilt.HYDRO_INFILT(hydroInfilt, Id_Select, KunsatModel_Infilt, N_SoilSelect)

			table.infilt.INFILT(Id_Select, N_SoilSelect, infiltOutput)
		end

	# PRINT OUTPUT ======================================================================================
	if option.globalopt.Ploting && !option.globalopt.Hypix
	println("		=== START: PLOTTING  ===")
	
		if option.smap.Plot_Kunsat  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			plotSmap.PLOT_KUNSAT(hydro, N_SoilSelect, smap; N_Se= 1000)
		end

		if option.globalopt.θΨ ≠ :No && option.hydro.Plot_θΨ # <>=<>=<>=<>=<>
			# plot.lab.σ_Ψm(hydro)
			if option.globalopt.Smap
				plotSmap.makie.HYDROPARAM(Ψ_θΨ, Ψ_KΨ, θ_θΨ, N_θΨ, N_SoilSelect, N_KΨ, K_KΨ, Id_Select, hydro, KunsatModel_Lab; smap=smap)

				# plotSmap.HYDROPARAM(Ψ_θΨ, Ψ_KΨ, θ_θΨ, N_θΨ, N_SoilSelect, N_KΨ, K_KΨ, Id_Select, hydro, KunsatModel_Lab; N_Se=1000, smap=[])
			else
				plot.lab.makie.HYDROPARAM(Ψ_θΨ, Ψ_KΨ, θ_θΨ, N_θΨ, N_SoilSelect, N_KΨ, K_KΨ, Id_Select, hydro, KunsatModel_Lab)
			end	
		end # option.globalopt.θΨ
		if  option.globalopt.θΨ ≠ :No && option.hydro.Plot_σ_Ψm && option.hydro.HydroModel == :Kosugi
			plot.lab.σ_Ψm(hydro)
		end
		if option.globalopt.Psd && option.psd.Plot_θr # <>=<>=<>=<>=<>
			plot.psd.PLOT_θr(∑Psd, N_SoilSelect, hydro, hydroPsd)
		end
		if option.globalopt.Psd && option.psd.Plot_IMP_Model # <>=<>=<>=<>=<>
			plot.psd.PLOT_IMP_MODEL(Id_Select, Rpart, N_Psd, ∑Psd, Psd, N_SoilSelect, hydro, paramPsd) 
		end
		if  option.globalopt.Psd && option.psd.Plot_Psd_θΨ && !option.psd.HydroParam
			println("			~ PSD WARNING Sorry cannot plot Plot_Psd_θΨ as option.psd.HydroParam==false ~")
		end
		if option.globalopt.Psd && option.psd.Plot_Psd_θΨ && option.psd.HydroParam # <>=<>=<>=<>=<>
			plot.psd.PLOT_PSD_θΨ(Ψ_θΨ, Ψ_Rpart, θ_θΨ, θ_Rpart, N_θΨ, N_SoilSelect, N_Psd, Id_Select, hydroPsd, hydro)
		end
		if option.globalopt.Infilt && option.infilt.Plot_∑Infilt  # <>=<>=<>=<>=<>
			plot.infilt.PLOT_∑INFILT(Id_Select, N_Infilt, N_SoilSelect, ∑Infilt_Obs, Tinfilt, ∑Infilt_3D, ∑Infilt_1D, infiltOutput)
		end
		# if option.globalopt.Infilt && option.infilt.Plot_SeIni_Range # <>=<>=<>=<>=<>
		# Removing GRUtils software to avoid conflict
		# 	# plot.infilt.PLOT_∑INFILT_SEINI(hydroInfilt, Id_Select, infiltOutput, infiltParam, N_SoilSelect)
		# end
		if  option.globalopt.Infilt && option.globalopt.θΨ ≠ :No && option.infilt.Plot_Sorptivity_SeIni # <>=<>=<>=<>=<>
			plot.infilt.PLOT_SORPTIVITY_SEINI(hydro, Id_Select, N_SoilSelect) 
		end
		if option.globalopt.Infilt && option.infilt.Plot_θΨ
			if option.globalopt.θΨ ≠ :No
				plot.infilt.PLOT_∑INFILT_θΨ(hydroInfilt, Id_Select, N_SoilSelect; hydro=hydro)
			else
				plot.infilt.PLOT_∑INFILT_θΨ(hydroInfilt, Id_Select, N_SoilSelect)
			end # option.globalopt.θΨ
		end # option.globalopt.Infilt

	println("=== END: PLOTTING  === \n")
	end # if option.globalopt.Ploting

	# Playing sounds...
		println("\007")

end  # function: START_TOOLBOX
# ..............................................................

println("\n\n===== START SOIL WATER TOOLBOX =====")
	@time START_TOOLBOX()
println("==== END SOIL WATER TOOLBOX ====")