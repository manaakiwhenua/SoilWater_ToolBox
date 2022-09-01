# =============================================================
#		module: tableSmap
# =============================================================
module tableSmap
   import ..path, ..tool, ..option, ..param, ..wrc
   import CSV, Tables, DataFrames, DelimitedFiles


   	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : θΨK
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         function θΨK(hydro, hydroOther, Id_Select, KunsatModel_Lab, N_SoilSelect, smap)
            println("    ~  $(path.Table_θΨK) ~")

            Matrix, FieldName_String = tool.readWrite.STRUCT_2_FIELDNAME(N_SoilSelect, hydro)

            Matrix2, FieldName_String2 = tool.readWrite.STRUCT_2_FIELDNAME(N_SoilSelect, hydroOther)

            # Concatenating matrices
               Matrix = hcat(Matrix, Matrix2)

               Matrix = hcat(Matrix, KunsatModel_Lab)

               FieldName_String = vcat(FieldName_String, FieldName_String2)

               pushfirst!(FieldName_String, string("Depth"))
               pushfirst!(FieldName_String, string("SoilName"))
               pushfirst!(FieldName_String, string("Id")) 
               push!(FieldName_String, string("KunsatModel"))
               
            CSV.write(path.Table_θΨK,Tables.table([Id_Select smap.Soilname[1:N_SoilSelect] smap.Depth[1:N_SoilSelect] Matrix]), header=FieldName_String )
      
         return nothing
         end  # function:  θΨK

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : Smap
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   """
   TopnetModel = ["ThetaS_Ch[mm3 mm-3]";"ThetaR_Ch[mm3 mm-3]";"LambdaCh_Ch[-]";"Hch_Ch[mm]";" Hga_Ch[mm]";"Ks_Vg[mm s-1]; "0mm"; "500mm"; "1000mm"; "2000mm"; "4000mm"; "10000mm"; "150000mm"]

   JulesModel_CH =  ["ThetaS_Ch[mm3 mm-3]";"ThetaR_Ch[mm3 mm-3]";"LambdaCh_Ch[-]";"Hch_Ch[mm]"; "Ks_Ch[mm s-1]";" Hga_Ch[mm]";"3300mm";"10000mm" ;"150000mm"]

   JulesModel_VangenuchtenJules = ["ThetaS_VgJules[mm3 mm-3]";"ThetaR_VgJules[mm3 mm-3]";"n_VgJules[-]";"Hvg_VgJules[mm]"; "Ks_VgJules[mm s-1]";"3300mm";"10000mm"]

   """
      function SMAP(Id_Select, N_SoilSelect, smap)
         println("    ~  $(path.Table_Smap) ~")

         # Header
            HeaderSmap = false # <true> the greek characters are replaced by alphabet; <false> original parameter names with no units usefull to use values in SoilWater-ToolBox

         # User input
            Option_BrooksCorey       = true
            Option_ClappHornberger   = true
            Option_VanGenuchten      = true
            Option_VanGenuchtenJules = true
            Option_Kosugi            = true
            Option_Kosugi_Table      = true

         Header = ["Id"; "SoilName"; "Depth_mm"; "IsTopsoil"; "RockFragment_%";"RockDepth_mm"; "MaxRootingDepth_mm"]

         Data = []
      
      # Select data
         # HydroModel_θΨ == "BrooksCorey"
         if Option_BrooksCorey # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
            HydroModel_θΨ = "BrooksCorey"

            Path_θΨ =  path.FileSoilHydro_Table₁ * string(HydroModel_θΨ) *  "_" * string(option.hydro.σ_2_Ψm) *  "_" * path.Model_Name * "_" * path.Table_θΨK₀
            
            if isfile(Path_θΨ)
               Select_θΨ = ["θs";"θr";"λbc";"Ψbc"; "Ks"; "Ψga"]
               
               Data_θΨ = Tables.matrix(CSV.File(Path_θΨ, select=Select_θΨ))

               try
                  Data = hcat(Data[1:N_SoilSelect, :], Data_θΨ[1:N_SoilSelect, :])
               catch
                  Data = Data_θΨ[1:N_SoilSelect, :]
               end
               
               if HeaderSmap
                  Header_θΨ = ["ThetaS_Bc[mm3 mm-3]";"ThetaR_Bc[mm3 mm-3]";"LambdaBc_Bc[-]";"Hbc_Bc[mm]";"Ks_Bc[mm s-1]";"Hga_Bc[mm]"]
               else
                  Header_θΨ = Select_θΨ
               end

               Header =  append!(Header, Header_θΨ)
            else
               @warn("\n \n WARNING Smap_Output: model simulation not found: $HydroModel_θΨ \n")
            end # if isfile(Path_θΨ)
         end # Option_BrooksCorey


         # HydroModel_θΨ == "ClappHornberger"
         if  Option_ClappHornberger # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
            HydroModel_θΨ = "ClappHornberger"

            Path_θΨ =  path.FileSoilHydro_Table₁ * string(HydroModel_θΨ) *  "_" * string(option.hydro.σ_2_Ψm) *  "_" * path.Model_Name * "_" * path.Table_θΨK₀ 

            if isfile(Path_θΨ)

               Select_θΨ = ["θs";"θr";"λch";"Ψch";"Ks";"Ψga"]

               Data_θΨ = Tables.matrix(CSV.File(Path_θΨ, select=Select_θΨ))

               try
                  Data = hcat(Data[1:N_SoilSelect, :], Data_θΨ[1:N_SoilSelect, :])
               catch
                  Data = Data_θΨ[1:N_SoilSelect, :]
               end

               if HeaderSmap
                  Header_θΨ = ["ThetaS_Ch[mm3 mm-3]";"ThetaR_Ch[mm3 mm-3]";"LambdaCh_Ch[-]";"Hch_Ch[mm]"; "Ks_Ch[mm s-1]";" Hga_Ch[mm]"]
               else
                  Header_θΨ = Select_θΨ
               end

               Header =  append!(Header, Header_θΨ)
            else
               @warn("\n \n WARNING Smap_Output: model simulation not found: $HydroModel_θΨ \n")
            end # if isfile(Path_θΨ)
         end #  Option_ClappHornberger

         
         # HydroModel_θΨ == "Vangenuchten"
         if Option_VanGenuchten # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
            HydroModel_θΨ = "Vangenuchten"

            Path_θΨ =  path.FileSoilHydro_Table₁ * string(HydroModel_θΨ) *  "_" * string(option.hydro.σ_2_Ψm) *  "_" * path.Model_Name * "_" * path.Table_θΨK₀ 

            if isfile(Path_θΨ)
               Select_θΨ = ["θs";"θr";"N";"Ψvg"; "Ks"]
               
               Data_θΨ = Tables.matrix(CSV.File(Path_θΨ, select=Select_θΨ))
                  
               try
                  Data = hcat(Data[1:N_SoilSelect, :], Data_θΨ[1:N_SoilSelect, :])
               catch
                  Data = Data_θΨ[1:N_SoilSelect, :]
               end
         
               if HeaderSmap
                  Header_θΨ = ["ThetaS_Vg[mm3 mm-3]";"ThetaR_Vg[mm3 mm-3]";"N_Vg[-]";"Hvg_Vg[mm]"; "Ks_Vg[mm s-1]"]
               else
                  Header_θΨ = Select_θΨ
               end

               Header =  append!(Header, Header_θΨ)
            else
               @warn("\n \n WARNING Smap_Output: model simulation not found: $HydroModel_θΨ \n")
            end # if isfile(Path_θΨ)
         end # Option_VanGenuchten


         # HydroModel_θΨ == "VangenuchtenJules"
         if Option_VanGenuchtenJules # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
            HydroModel_θΨ = "VangenuchtenJules"

            Path_θΨ =  path.FileSoilHydro_Table₁ * string(HydroModel_θΨ) *  "_" * string(option.hydro.σ_2_Ψm) *  "_" * path.Model_Name * "_" * path.Table_θΨK₀ 

            if isfile(Path_θΨ)
               Select_θΨ = ["θs";"θr";"N";"Ψvg"; "Ks"]
               
               Data_θΨ = Tables.matrix(CSV.File(Path_θΨ, select=Select_θΨ))
                  
               try
                  Data = hcat(Data[1:N_SoilSelect, :], Data_θΨ[1:N_SoilSelect, :])
               catch
                  Data = Data_θΨ[1:N_SoilSelect, :]
               end
         
               if HeaderSmap
                  Header_θΨ = ["ThetaS_VgJules[mm3 mm-3]";"ThetaR_VgJules[mm3 mm-3]";"n_VgJules[-]";"Hvg_VgJules[mm]"; "Ks_VgJules[mm s-1]"]
               else
                  Header_θΨ = Select_θΨ
               end

               Header =  append!(Header, Header_θΨ)
            else
               @warn("\n \n WARNING Smap_Output: model simulation not found: $HydroModel_θΨ \n")
            end # if isfile(Path_θΨ)
         end # Option_VanGenuchten


         # HydroModel_θΨ == "Kosugi"
         if Option_Kosugi # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
            HydroModel_θΨ = "Kosugi"

            Path_θΨ =  path.FileSoilHydro_Table₁ * string(HydroModel_θΨ) *  "_" * string(option.hydro.σ_2_Ψm) *  "_" * path.Model_Name * "_" * path.Table_θΨK₀ 

            if isfile(Path_θΨ)
               Select_θΨ =["θs";"θr";"Ks";"Ψm";"σ";"σMac";"ΨmMac";"θsMacMat"; "θsMacMat_ƞ"]
                   
               Data_θΨ = Tables.matrix(CSV.File(Path_θΨ, select=Select_θΨ))
                  
               try
                  Data = hcat(Data[1:N_SoilSelect, :], Data_θΨ[1:N_SoilSelect, :])
               catch
                  Data = Data_θΨ[1:N_SoilSelect, :]
               end
         
               if HeaderSmap
                  Header_θΨ = ["ThetaS_Kg[mm3 mm-3]";"ThetaR_Kg[mm3 mm-3]";"Ks_Kg[mm s-1]";"Hm_Kg[mm]";"Sigma_Kg";"SigmaMac_Kg";"HmMac_Kg[mm]";"ThetaSMacMat_Kg[mm3 mm-3]"]
               else
                  Header_θΨ = Select_θΨ
               end

               Header =  append!(Header, Header_θΨ)
            else
               @warn("\n \n WARNING Smap_Output: model simulation not found: $HydroModel_θΨ \n")
            end # if isfile(Path_θΨ)
         end # Option_Kosugi


      # HydroModel_θΨ == "Option_Kosugi_Table"
      if Option_Kosugi_Table # <>=<>=<>=<>=<>=<>=<>=<>=<>=<>
         # HydroModel_θΨ = "Kosugi"

         Path_θΨ =  path.Table_KosugiθΨ

         if isfile(Path_θΨ)
            Select_θΨ = string.(Int64.(param.hydro.smap.Ψ_Table)) .* "mm"
            
            Data_θΨ = Tables.matrix(CSV.File(Path_θΨ, select=Select_θΨ))
               
            try
               Data = hcat(Data[1:N_SoilSelect, :], Data_θΨ[1:N_SoilSelect, :])
            catch
               Data = Data_θΨ[1:N_SoilSelect, :]
            end
      
            Header_θΨ = Select_θΨ

            Header =  append!(Header, Header_θΨ)
         else
            @warn("\n \n WARNING Smap_Output: model simulation not found: $HydroModel_θΨ \n")
         end # if isfile(Path_θΨ)
      end # Option_Kosugi

         
      # COMBINING OUTPUTS   
         Output = Tables.table([string.(Id_Select[1:N_SoilSelect]) smap.Soilname[1:N_SoilSelect] smap.Depth[1:N_SoilSelect] smap.IsTopsoil[1:N_SoilSelect] smap.RockFragment[1:N_SoilSelect] smap.RockDepth[1:N_SoilSelect] smap.MaxRootingDepth[1:N_SoilSelect] Data[1:N_SoilSelect,:]])
      
         CSV.write(path.Table_Smap, Output, header=Header)	

      return nothing
      end  # function:  smap
	# ............................................................

  
end  # module: tableSmap
# ............................................................