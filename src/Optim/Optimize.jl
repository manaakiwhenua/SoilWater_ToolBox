 # =============================================================
#		MODULE: optimize
# =============================================================
 module optimize

   import ..option
   export SEARCHRANGE
   
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : SEARCHRANGE
	#		Required by BlackBoxOptim
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function SEARCHRANGE(optim)
         ParamOpt_Min₂ = copy(optim.ParamOpt_Min)
         ParamOpt_Max₂ = copy(optim.ParamOpt_Max)

         # Log transform
         for iParam=1:optim.NparamOpt
            if optim.ParamOpt_LogTransform[iParam]
               ParamOpt_Min₂[iParam] = log1p(optim.ParamOpt_Min[iParam])
               ParamOpt_Max₂[iParam] = log1p(optim.ParamOpt_Max[iParam])
            end
         end

      # Making sure that for constrained optimisation Ψm is between 0 & 1
         if (option.hydro.σ_2_Ψm==:Constrained) && ("Ψm" ∈ optim.ParamOpt)
            iψm = findfirst(isequal("Ψm"), optim.ParamOpt)[1]

            ParamOpt_Min₂[iψm] = 0.0
            ParamOpt_Max₂[iψm] = 1.0
         end

         SearchRange = (collect(zip(Float64.(ParamOpt_Min₂), Float64.(ParamOpt_Max₂))))

      return SearchRange
      end  # function: SEARCHRANGE
end # module optimize