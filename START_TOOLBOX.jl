##========================================================================================
##                                                                                      ##
##                                 Soil Water ToolBox                                   ##
##                                                                                     ##
##========================================================================================

using Suppressor
@suppress_err begin
	include("Option.jl")
	# Install packages to run program
		if option.DownloadPackage
			include("Packages.jl")
		end
	include("Path.jl")
	include("Cst.jl")
	include("Param.jl")
	include("Tool.jl")
	include("Read.jl")
	include("Hydro\\TotalPorosity.jl")
	include("Hydro\\HydroStruct.jl")
	include("Hydro\\HydroInitialize.jl")
	include("Hydro\\WaterRetentionCurve.jl")
	include("Hydro\\Kunsat.jl")
	include("Stats.jl")
	
	include("Hydro\\ObjectiveFunction_Hydro.jl")
	include("Hydro\\START_Hydro.jl")
	include("Hydro\\HydroRelation.jl")
	include("Psd\\PsdThetar.jl")
		if option.Infilt
			if option.infilt.SortivityVersion == "NonInfinity"
				include("Infilt\\SorptivityNonInfinity.jl")
			elseif option.infilt.SortivityVersion == "Traditional"
				include("Infilt\\SorptivityTraditional.jl")
			end
			include("Infilt\\Best_Univ.jl")
			include("Infilt\\QuasiExact.jl")
			include("Infilt\\TimeTransSteady.jl")
			include("Infilt\\InfiltStruct.jl")
			include("Infilt\\InfiltInitialize.jl")
			include("Infilt\\START_Infilt.jl")
			# include("Infilt\\OptInfilt.jl")
		end
		if option.Psd
			include("Psd\\PsdStruct.jl")
			include("Psd\\PsdInitialize.jl")
		end
		if option.Psd
			include("Psd\\PsdFunc.jl")
			include("Psd\\PsdOpt.jl")
			include("Psd\\START_PSD.jl")
		end
		if option.Plot
			include("Plot.jl")
		end
	include("Table.jl")
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		FUNCTION : START_TOOLBOX
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function START_TOOLBOX()

	println("=== START: READING ===")
		# Selecting soils of interest 
		Id_Select, N_SoilSelect = read.ID()

		# Reinforcing the maximum of iSoil to simulate
			N_SoilSelect = Int(min(N_SoilSelect, param.N_iSoil_Simulations))
			Id_Select = Id_Select[1:N_SoilSelect]

		if option.???? ??? "No"
			??_????, ??_????, N_???? = read.????(Id_Select, N_SoilSelect)

			RockW_????, ??_Rock_????, ??bSoil_????, ??p_Fine_???? = read.??_????(Id_Select, N_SoilSelect)
		end # option.???? ??? "No"

		if option.hydro.Kunsat??
			K_K??, ??_K??, N_K?? = read.KUNSAT??(Id_Select, N_SoilSelect)
		end # option.hydro.Kunsat??

		if option.Psd
			Rpart, ???Psd, N_Psd = read.PSD(Id_Select, N_SoilSelect)

			RockW_Psd, ??_Rock_Psd, ??bSoil_Psd, ??p_Fine_Psd = read.??_PSD(Id_Select, N_SoilSelect)
		else
			???Psd = zeros(Float64, N_SoilSelect, 1)
			Rpart = zeros(Float64, N_SoilSelect, 1)
			N_Psd= zeros(Float64, N_SoilSelect)
		end # option.Psd
		
		if option.Infilt
			Tinfilt, ???Infilt_Obs, N_Infilt, infiltParam  = read.INFILTRATION(Id_Select, N_SoilSelect)

			RockW_Infilt, ??_Rock_Infilt, ??bSoil_Infilt, ??p_Fine_Infilt = read.??_INFILTRATION(Id_Select, N_SoilSelect)
		end


	println("=== END  : READING === \n")


	if option.???? ??? "No"
	println("=== START: DERIVING HYDRO PARAMETERS  ===")
		# INITIALIZES HYDRAULIC PARAMETERS STRUCT INDEPENDENTLY OF THE SELECTED MODEL
			hydro = hydroStruct.HYDROSTRUCT(N_SoilSelect)

		# Total Porosity= ??
			hydro.?? = ??.??B_2_??(N_SoilSelect, RockW_????, ??_Rock_????, ??bSoil_????, ??p_Fine_????)

		if option.hydro.Kunsat??
			# Structure of hydro
			hydro = hydroParam.START_HYDROPARAM(N_SoilSelect=N_SoilSelect, ???Psd=???Psd, ??_????=??_????, ??_????=??_????, N_????=N_????, K_K??=K_K??, ??_K??=??_K??, N_K??=N_K??, hydro=hydro, optionHydro=option.hydro)
		else
			hydro = hydroParam.START_HYDROPARAM(N_SoilSelect=N_SoilSelect, ???Psd=???Psd, ??_????=??_????, ??_????=??_????, N_????=N_????, hydro=hydro, optionHydro=option.hydro)
		end
	println("=== END  : DERIVING HYDRO PARAMETERS  === \n")
	else
		hydro = []
	end


	if option.Psd
	println("=== START: PSD MODEL  ===")
		# Structure of hydroPsd
			hydroPsd = hydroStruct.HYDROSTRUCT(N_SoilSelect)

		# Total Porosity= ??
			hydroPsd.?? = ??.??B_2_??(N_SoilSelect,RockW_Psd, ??_Rock_Psd, ??bSoil_Psd, ??p_Fine_Psd)

			paramPsd, N_Psd, ??_Rpart, ??_Rpart, Psd, hydroPsd = psd.START_PSD(N_SoilSelect, ??_????, ??_????, N_????, Rpart, ???Psd, N_Psd, hydro, hydroPsd)

		if  option.psd.HydroParam
			hydroPsd = hydroParam.START_HYDROPARAM(N_SoilSelect=N_SoilSelect, ???Psd=???Psd, ??_????=??_Rpart, ??_????=??_Rpart, N_????=N_Psd, hydro=hydroPsd, optionHydro=option.psd)
		end

	println("=== END : PSD MODEL  === \n")
	else
		??_Rpart = zeros(Float64, N_SoilSelect,1)
		??_Rpart = zeros(Float64, N_SoilSelect,1)
		hydroPsd = hydroStruct.HYDROSTRUCT(N_SoilSelect)
		N_Psd = zeros(Float64, N_SoilSelect)
	end

	
	if option.Infilt
	println("=== START: INFILTRATION  ===")
		# Structure of hydroInfilt
			hydroInfilt = hydroStruct.HYDROSTRUCT(N_SoilSelect)

		# Total Porosity= ??
			hydroInfilt.?? = ??.??B_2_??(N_SoilSelect, RockW_Infilt, ??_Rock_Infilt, ??bSoil_Infilt, ??p_Fine_Infilt)

		infiltOutput, hydroInfilt, ???Infilt = infilt.START_INFILTRATION(???Infilt_Obs, ???Psd, hydro, hydroInfilt, infiltParam, N_Infilt, N_SoilSelect, Tinfilt, Id_Select)

	println("=== END  : INFILTRATION  === \n")
	else
		hydroInfilt = []
	end

	println("=== START: WRITING TABLE  ===")
		if option.???? ??? "No"
			table.hydroParam.????K(Id_Select[1:N_SoilSelect], N_SoilSelect, hydro)
		end

		if option.Psd
			table.psd.PSD(Id_Select[1:N_SoilSelect], N_SoilSelect, paramPsd)

			if option.psd.HydroParam
				table.psd.????K_PSD(Id_Select, N_SoilSelect, hydroPsd)
			end	
		end

		if option.Infilt
			table.infilt.HYDRO_INFILT(Id_Select, N_SoilSelect, hydroInfilt)

			table.infilt.INFILT(Id_Select, N_SoilSelect, infiltOutput)
		end
	println("=== END  : WRITING TABLE  === \n")


	if option.Plot
	println("=== START: PLOTTING  ===")

		if option.???? ??? "No" && option.hydro.Plot_????
			plot.HYDROPARAM(??_????, ??_K??, ??_????, N_????, N_SoilSelect, N_K??, K_K??, Id_Select, hydro)
		end # option.Plot_WaterRetentionCurve

		if option.Psd && option.psd.Plot_??r
			plot.PLOT_??r(???Psd, N_SoilSelect, hydro, paramPsd)
		end

		if option.Psd && option.psd.Plot_IMP_model
			plot.PLOT_IMP_model(Id_Select, Rpart, N_Psd, ???Psd, Psd, N_SoilSelect, hydro, paramPsd) 
		end

		if option.infilt.Plot_???Infilt && option.Infilt
			plot.PLOT_???INFILT(Id_Select, N_Infilt, N_SoilSelect, ???Infilt_Obs, Tinfilt, ???Infilt, infiltOutput)

			# plot.PLOT_TREANSSTEADY(Id_Select, N_Infilt, N_SoilSelect, ???Infilt_Obs, Tinfilt, ???Infilt, infiltOutput)
		end

	println("=== END: PLOTTING  === \n")
	end # if option.Plot
		
end  # function: START_TOOLBOX

println("\n\n===== START SOIL WATER TOOLBOX ==== \n")
	@time START_TOOLBOX()
println("==== END SOIL WATER TOOLBOX ===")