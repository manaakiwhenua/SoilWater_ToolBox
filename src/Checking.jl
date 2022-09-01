# =============================================================
#		module: checking
# =============================================================
module checking
   import ..path, ..option

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : CHECKING
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function CHECKING(optim)

         	# CHECKING FOR UNCONSISTENCY WITH OPTIONS
			if option.hydro.HydroModel==:Kosugi && "θsMacMat" ∈ optim.ParamOpt
				error("*** option.hydro.HydroModel==:Kosugi && option.hydro.HydroModel==:Bimodal THAN option.hydro.HydroModel == :Φ ***")
					
			elseif option.hydro.θrOpt≠:Opt && "θr" ∈ optim.ParamOpt
				error("*** option.hydro.θrOpt≠:Opt && θr ∈ ParamOpt THAN do not optimize θr ***")
			              
			elseif option.hydro.σ_2_Ψm ==:UniqueRelationship && "Ψm" ∈ optim.ParamOpt
				error("*** option.hydro.σ_2_Ψm ==:UniqueRelationship THAN Ψm does not need to be optimised ***")
			
			elseif option.hydro.HydroModel == :Kosugi && (option.hydro.σ_2_Ψm ==:Constrained && "Ψm" ∉ optim.ParamOpt) 
				error("*** option.hydro.σ_2_Ψm ==:Constrained THAN Ψm needs to be optimised ***")

			elseif option.hydro.KsOpt == :Data && "Ks" ∈ optim.ParamOpt
				error("*** option.hydro.KsOpt == :Data THAN Ks does not need to be optimized ***")

			elseif  (option.hydro.θrOpt==:σ_2_θr) && ("θr" ∈ optim.ParamOpt)
				error("*** option.hydro.θrOpt==:σ_2_θr THAN θr does not need to be optimized ***") 

			elseif (option.hydro.θrOpt==:σ_2_θr) && ("σ" ∉ optim.ParamOpt)
				error("*** option.hydro.θrOpt==:σ_2_θr THAN σ needs to be optimized ***")
				
			elseif  (option.hydro.θrOpt==:ParamPsd) && ("θr"∈ optim.ParamOpt) # Derive θr frpm PSD
				error("*** option.Hydro.θrOpt==:ParamPsd && θr does not need to be optimized ***")

			elseif  (option.hydro.θrOpt==:ParamPsd) && ("θr"∈ optim.ParamOpt) # Derive θr frpm PSD
				error("*** option.Hydro.θrOpt==:ParamPsd && θr does not need to be optimized ***")

			elseif  (option.hydro.θrOpt==:ParamPsd) && ("θr"∉ optim.ParamOpt) && !(option.globalopt.Psd) # Derive θr frpm PSD
				error("*** option.Hydro.θrOpt==:ParamPsd THAN option.globalopt.Psd=true ***")
			
			elseif ("Ks" ∈ optim.ParamOpt) && (option.hydro.KunsatΨ==false)
				error("***  (Ks ∈ optim.ParamOpt) && (KunsatΨ==false) THAN option.KunsatΨ=true ***")

         elseif option.smap.UsePointKosugiBimodal && option.hydro.KunsatΨ && "Ks" ∉ optim.ParamOpt
            error("***  Ks  ∉ optim.ParamOpt && option.smap.UsePointKosugiBimodal && option.hydro.KunsatΨ THAN Ks ∈ optim.ParamOpt***")
			
			elseif  option.hydro.HydroModel==:Kosugi && (option.smap.AddPointKosugiBimodal) && option.hydro.KunsatΨ && option.globalopt.Smap
				error("option.hydro.HydroModel==:Kosugi && (option.smap.AddPointKosugiBimodal) THEN option.hydro.KunsatΨ=false OR UsePointKosugiBimodal = true")

			elseif  option.hydro.HydroModel==:Kosugi && (option.smap.AddPointKosugiBimodal) &&  "Ks" ∈ optim.ParamOpt && option.Smap
				error("option.hydro.HydroModel==:Kosugi && (option.smap.AddPointKosugiBimodal) THEN Ks ∉ optim.ParamOpt")
			end # Check error

         return nothing
      end  # function: CHECKING
   
end  # module: checking
# ............................................................