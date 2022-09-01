# =============================================================
#		module: stoneCorrection
# =============================================================
module stoneCorrection
   using Polynomials
   export STONECORRECTION

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : STONECORRECTION_NONEWETABLE
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function STONECORRECTION(N_SoilSelect, N_θΨ, smap, θ_θΨ, Ψ_θΨ)

         for iZ = 1:N_SoilSelect 
            for iθ=1:N_θΨ[iZ]
               θ_θΨ[iZ,iθ] =  θ_θΨ[iZ,iθ] * (1.0 - smap.RockFragment[iZ])
            end # for iθ=1:N_θΨ[iZ]
         end #  for iZ = 1:N_SoilSelect
         
      return  θ_θΨ
      end  # function: STONECORRECTION_NONEWETABLE


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : STONECORRECTION_NONEWETABLE_HYDRO
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function STONECORRECTION_HYDRO(hydroParam, N_SoilSelect, optim, smap)

         for iZ = 1:N_SoilSelect
            
            hydroParam.θs_Max[iZ] =  hydroParam.θs_Max[iZ] * (1.0 - smap.RockFragment[iZ])
            hydroParam.θs_Min[iZ] =  hydroParam.θs_Min[iZ] * (1.0 - smap.RockFragment[iZ])
            
            if ("θs" ∈ optim.ParamOpt)
               iθs = findfirst(isequal("θs"), optim.ParamOpt)[1]
               optim.ParamOpt_Min[iθs] = hydroParam.θs_Min[iZ]
               optim.ParamOpt_Max[iθs] = hydroParam.θs_Max[iZ]
            end # "θs"

            hydroParam.θr_Max[iZ] =  hydroParam.θr_Max[iZ] * (1.0 - smap.RockFragment[iZ])
            hydroParam.θr_Min[iZ] =  hydroParam.θr_Min[iZ] * (1.0 - smap.RockFragment[iZ])
            
            if ("θr" ∈ optim.ParamOpt)
               iθr = findfirst(isequal("θr"), optim.ParamOpt)[1]
               optim.ParamOpt_Min[iθr] = hydroParam.θr_Min[iZ]
               optim.ParamOpt_Max[iθr] = hydroParam.θr_Max[iZ]
            end # "θr"

            hydroParam.Ks_Max[iZ] =   hydroParam.Ks_Max[iZ] * (1.0 - smap.RockFragment[iZ])
            hydroParam.Ks_Min[iZ] =   hydroParam.Ks_Min[iZ] * (1.0 - smap.RockFragment[iZ])
            
            if ("Ks" ∈ optim.ParamOpt)
               iKs = findfirst(isequal("Ks"), optim.ParamOpt)[1]
               optim.ParamOpt_Min[iKs] = hydroParam.Ks_Min[iZ]
               optim.ParamOpt_Max[iKs] = hydroParam.Ks_Max[iZ]
            end # "Ks"
         end #  for iZ = 1:N_SoilSelect
         
      return  hydroParam, optim
      end  # function: STONECORRECTION_NONEWETABLE



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		FUNCTION : STONECORRECTION_WETABLE
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   function STONECORRECTION_WETTABLE(N_SoilSelect, N_θΨ, rfWetable, smap, θ_θΨ, Ψ_θΨ)
      for iZ = 1:N_SoilSelect
         
         iRockClass = rfWetable.RockClass_Dict[smap.RockClass[iZ]]

         WettableFunction = Polynomial(rfWetable.RockClass_Polynomial_Array[iRockClass][1])

         for iθ=1:N_θΨ[iZ]
            θ_Wettable = WettableFunction(log1p(Ψ_θΨ[iZ,iθ]))

            θ_θΨ[iZ,iθ] =  θ_θΨ[iZ,iθ] + θ_Wettable * smap.RockFragment[iZ]
         end # for iθ=1:N_θΨ[iZ]
      end #  for iZ = 1:N_SoilSelect
      
   return  θ_θΨ
   end  # function: STONECORRECTION_NONEWETABLE

end  # module: stoneCorrection
# ............................................................