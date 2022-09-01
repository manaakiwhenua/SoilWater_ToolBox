# =============================================================
#		MODULE: infiltStruct
# =============================================================
module infiltStruct
	import ..option, ..tool
	export INFILTSTRUCT

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		STRUCTURE : HYDRAULIC
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		mutable struct INFILT
            Sorptivity         :: 	Vector{Float64}
            iT_TransSteady_Data :: 	Vector{Int64}
            T_TransSteady_Data  :: 	Vector{Float64}
            Nse_Trans          ::	Vector{Float64}
				Nse_Steady         ::	Vector{Float64}
				Nse			       ::	Vector{Float64}
		end # struct KOSUGI

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INFILTSTRUCT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INFILTSTRUCT(N_SoilSelect)
         FieldName           = Array{Symbol}(undef, 1) # Need to put
         Sorptivity          = zeros(Float64, N_SoilSelect)
         iT_TransSteady_Data = zeros(Int64, N_SoilSelect)
         T_TransSteady_Data  = zeros(Float64, N_SoilSelect)
         Nse_Trans           = zeros(Float64, N_SoilSelect)
         Nse_Steady          = zeros(Float64, N_SoilSelect)
         Nse                 = zeros(Float64, N_SoilSelect)
			
			return infiltOutput = INFILT(Sorptivity, iT_TransSteady_Data, T_TransSteady_Data, Nse_Trans, Nse_Steady, Nse)

		end  # function: INFILTSTRUCT

end # module infiltStruct
# ............................................................