# =============================================================
#		MODULE: hydroStruct
# =============================================================
module hydroStruct
	import ..option, ..tool
	export HYDROSTRUCT

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		STRUCTURE : HYDRAULIC
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		mutable struct KOSUGI  # <>=<>=<>=<>=<>=<>=<>=<>=<>
         θs             :: Vector{Float64}
         θr             :: Vector{Float64}
         Ks             :: Vector{Float64}
         Ψm             :: Vector{Float64}
         σ              :: Vector{Float64}
         θsMacMat_ƞ     ::	Vector{Float64}
         σMac           :: Vector{Float64}
         ΨmMac          ::	Vector{Float64}
         So             :: Vector{Float64}
         θsMacMat       ::	Vector{Float64}
         Φ              ::	Vector{Float64}
			
         θs_Min         :: Vector{Float64}
         θr_Min         :: Vector{Float64}
         Ks_Min         :: Vector{Float64}
         Ψm_Min         :: Vector{Float64}
         σ_Min          :: Vector{Float64}
         θsMacMat_ƞ_Min ::	Vector{Float64}
         σMac_Min       :: Vector{Float64}
         ΨmMac_Min      ::	Vector{Float64}
         So_Min         :: Vector{Float64}
         θsMacMat_Min   ::	Vector{Float64}
         Φ_Min          ::	Vector{Float64}

         θs_Max         :: Vector{Float64}
         θr_Max         :: Vector{Float64}
         Ks_Max         :: Vector{Float64}
         Ψm_Max         :: Vector{Float64}
         σ_Max          :: Vector{Float64}
         θsMacMat_ƞ_Max ::	Vector{Float64}
         σMac_Max       :: Vector{Float64}
         ΨmMac_Max      ::	Vector{Float64}
         So_Max         :: Vector{Float64}
         θsMacMat_Max   ::	Vector{Float64}
         Φ_Max          ::	Vector{Float64}
		end # struct KOSUGI

		
		mutable struct VANGENUCHTEN # <>=<>=<>=<>=<>=<>=<>=<>=<>
         θs      ::	Vector{Float64}
         θr      ::	Vector{Float64}
         N       ::	Vector{Float64}
         Ψvg     ::	Vector{Float64}
         Ks      ::	Vector{Float64}
         Km      ::	Vector{Float64}
         Φ       ::	Vector{Float64}
			
         θs_Min  ::	Vector{Float64}
         θr_Min  ::	Vector{Float64}
         N_Min   ::	Vector{Float64}
         Ψvg_Min ::	Vector{Float64}
         Ks_Min  ::	Vector{Float64}
         Km_Min  ::	Vector{Float64}
         Φ_Min   ::	Vector{Float64}
			
         θs_Max  ::	Vector{Float64}
         θr_Max  ::	Vector{Float64}
         N_Max   ::	Vector{Float64}
         Ψvg_Max ::	Vector{Float64}
         Ks_Max  ::	Vector{Float64}
         Km_Max  ::	Vector{Float64}
         Φ_Max   ::	Vector{Float64}
		end # struct VANGENUCHTEN


		mutable struct BROOKS_COREY # <>=<>=<>=<>=<>=<>=<>=<>=<>
         θs      ::	Vector{Float64}
         θr      ::	Vector{Float64}
         λbc     ::	Vector{Float64}
         Ψbc     ::	Vector{Float64}
         Ks      ::	Vector{Float64}
         Φ       ::	Vector{Float64}
         Ψga     ::	Vector{Float64}
			
         θs_Min  ::	Vector{Float64}
         θr_Min  ::	Vector{Float64}
         λbc_Min ::	Vector{Float64}
         Ψbc_Min ::	Vector{Float64}
         Ks_Min  ::	Vector{Float64}
         Φ_Min   ::	Vector{Float64}
			
         θs_Max  ::	Vector{Float64}
         θr_Max  ::	Vector{Float64}
         λbc_Max ::	Vector{Float64}
         Ψbc_Max ::	Vector{Float64}
         Ks_Max  ::	Vector{Float64}
         Φ_Max   ::	Vector{Float64}
		end # struct BROOKS COREY


		mutable struct CLAPP_HORNBERGER # <>=<>=<>=<>=<>=<>=<>=<>=<>
			θs        ::	Vector{Float64}
         θr        ::	Vector{Float64}
         λch       ::	Vector{Float64}
         Ψch       ::	Vector{Float64}
         Ks        ::	Vector{Float64}
			Φ         ::	Vector{Float64}
         Ψga       ::	Vector{Float64}
			
			θs_Min        ::	Vector{Float64}
         θr_Min        ::	Vector{Float64}
         λch_Min       ::	Vector{Float64}
         Ψch_Min       ::	Vector{Float64}
         Ks_Min        ::	Vector{Float64}
			Φ_Min         ::	Vector{Float64}
			
			θs_Max        ::	Vector{Float64}
         θr_Max        ::	Vector{Float64}
         λch_Max       ::	Vector{Float64}
         Ψch_Max       ::	Vector{Float64}
         Ks_Max        ::	Vector{Float64}
         Φ_Max         ::	Vector{Float64}
		end # struct CLAPP_HORNBERGER


		mutable struct HYDRO_OTHER # <>=<>=<>=<>=<>=<>=<>=<>=<>
         Nse     :: 	Vector{Float64}
         Nse_θΨ  :: 	Vector{Float64}
         Nse_KΨ  :: 	Vector{Float64}
         Rmse    :: 	Vector{Float64}
         Rmse_θΨ :: 	Vector{Float64}
         Rmse_KΨ :: 	Vector{Float64}
		end


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : HYDROSTRUCT
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function HYDROSTRUCT(N_SoilSelect::Int64) # <>=<>=<>=<>=<>=<>=<>=<>=<>
			# For all models
         θs         = fill(0.0::Float64, N_SoilSelect)
         θr         = fill(0.0::Float64, N_SoilSelect)
         Ks         = fill(0.0::Float64, N_SoilSelect)
         θsMacMat   = fill(0.0::Float64, N_SoilSelect)
         θsMacMat_ƞ = fill(0.0::Float64, N_SoilSelect)
         Φ          = fill(0.0::Float64, N_SoilSelect)
         So         = fill(0.0::Float64, N_SoilSelect)

         θs_Min     = fill(0.0::Float64, N_SoilSelect)
         θr_Min     = fill(0.0::Float64, N_SoilSelect)
			Ks_Min     = fill(0.0::Float64, N_SoilSelect)
			Φ_Min      = fill(0.0::Float64, N_SoilSelect)
			So_Min     = fill(0.0::Float64, N_SoilSelect)
			
         θs_Max     = fill(0.0::Float64, N_SoilSelect)
         θr_Max     = fill(0.0::Float64, N_SoilSelect)
			Ks_Max     = fill(0.0::Float64, N_SoilSelect)
			Φ_Max      = fill(0.0::Float64, N_SoilSelect)
			So_Max     = fill(0.0::Float64, N_SoilSelect)
			
			if option.hydro.HydroModel == :Kosugi # <>=<>=<>=<>=<>
				σ              = fill(0.0::Float64, N_SoilSelect)
            Ψm             = fill(0.0::Float64, N_SoilSelect)
            σMac           = fill(0.0::Float64, N_SoilSelect)
            ΨmMac          = fill(0.0::Float64, N_SoilSelect)
				
            Ψm_Min         = fill(0.0::Float64, N_SoilSelect)
            σ_Min          = fill(0.0::Float64, N_SoilSelect)
            θsMacMat_ƞ_Min = fill(0.0::Float64, N_SoilSelect)
            σMac_Min       = fill(0.0::Float64, N_SoilSelect)
            ΨmMac_Min      = fill(0.0::Float64, N_SoilSelect)
            θsMacMat_Min   = fill(0.0::Float64, N_SoilSelect)

            Ψm_Max         = fill(0.0::Float64, N_SoilSelect)
            σ_Max          = fill(0.0::Float64, N_SoilSelect)
            θsMacMat_ƞ_Max = fill(0.0::Float64, N_SoilSelect)
            σMac_Max       = fill(0.0::Float64, N_SoilSelect)
            ΨmMac_Max      = fill(0.0::Float64, N_SoilSelect)
            θsMacMat_Max   = fill(0.0::Float64, N_SoilSelect)
          
				hydro = KOSUGI(θs, θr, Ks, σ, Ψm, θsMacMat_ƞ, σMac, ΨmMac, So, θsMacMat, Φ, θs_Min, θr_Min, Ks_Min, σ_Min,Ψm_Min, θsMacMat_ƞ_Min, σMac_Min, ΨmMac_Min, So_Min, θsMacMat_Min, Φ_Min, θs_Max, θr_Max, Ks_Max, σ_Max, Ψm_Max, θsMacMat_ƞ_Max, σMac_Max, ΨmMac_Max, So_Max, θsMacMat_Max, Φ_Max)
				return hydro

			elseif option.hydro.HydroModel == :Vangenuchten || option.hydro.HydroModel == :VangenuchtenJules # <>=<>=<>=<>=<>
            N       = fill(0.0::Float64, N_SoilSelect)
            Ψvg     = fill(0.0::Float64, N_SoilSelect)
            Km      = fill(0.0::Float64, N_SoilSelect)

            θs_Min  = fill(0.0::Float64, N_SoilSelect)
            θr_Min  = fill(0.0::Float64, N_SoilSelect)
            N_Min   = fill(0.0::Float64, N_SoilSelect)
            Ψvg_Min = fill(0.0::Float64, N_SoilSelect)
            Ks_Min  = fill(0.0::Float64, N_SoilSelect)
            Km_Min  = fill(0.0::Float64, N_SoilSelect)
			
            θs_Max  = fill(0.0::Float64, N_SoilSelect)
            θr_Max  = fill(0.0::Float64, N_SoilSelect)
            N_Max   = fill(0.0::Float64, N_SoilSelect)
            Ψvg_Max = fill(0.0::Float64, N_SoilSelect)
            Ks_Max  = fill(0.0::Float64, N_SoilSelect)
            Km_Max  = fill(0.0::Float64, N_SoilSelect)

				hydro = VANGENUCHTEN(θs, θr, N, Ψvg, Ks, Km, Φ, θs_Min,θr_Min, N_Min, Ψvg_Min, Ks_Min, Km_Min, Φ_Min, θs_Max, θr_Max, N_Max, Ψvg_Max, Ks_Max, Km_Max, Φ_Max) 
				return hydro


			elseif option.hydro.HydroModel == :BrooksCorey # <>=<>=<>=<>=<>=<>
            λbc     = fill(0.0::Float64, N_SoilSelect)
            Ψbc     = fill(0.0::Float64, N_SoilSelect)
            λbc_Min = fill(0.0::Float64, N_SoilSelect)
            Ψbc_Min = fill(0.0::Float64, N_SoilSelect)
            λbc_Max = fill(0.0::Float64, N_SoilSelect)
            Ψbc_Max = fill(0.0::Float64, N_SoilSelect)
            Ψga     = fill(0.0::Float64, N_SoilSelect)

				hydro = BROOKS_COREY(θs, θr, λbc, Ψbc, Ks, Φ, Ψga, θs_Min, θr_Min, λbc_Min, Ψbc_Min, Ks_Min, Φ_Min, θs_Max, θr_Max, λbc_Max, Ψbc_Max, Ks_Max, Φ_Max)
				return hydro

			elseif option.hydro.HydroModel == :ClappHornberger # <>=<>=<>=<>=<>=<>
				λch = fill(0.0::Float64, N_SoilSelect)
				Ψch  = fill(0.0::Float64, N_SoilSelect)
            Ψga  = fill(0.0::Float64, N_SoilSelect)

				λch_Min = fill(0.0::Float64, N_SoilSelect)
				Ψch_Min  = fill(0.0::Float64, N_SoilSelect)

				λch_Max = fill(0.0::Float64, N_SoilSelect)
				Ψch_Max  = fill(0.0::Float64, N_SoilSelect)

				hydro = CLAPP_HORNBERGER(θs, θr, λch, Ψch, Ks, Φ, Ψga,  θs_Min, θr_Min, λch_Min, Ψch_Min, Ks_Min, Φ_Max, θs_Max, θr_Max, λch_Max, Ψch_Max, Ks_Max, Φ_Max)
				return hydro
			end # option.hydro.HydroModel

		end #  function HYDROSTRUCT


		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#		FUNCTION : HYDRO_OTHER
		# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			function HYDRO_OTHERS(N_SoilSelect::Int64)
				Nse = fill(0.0::Float64, N_SoilSelect)
            Nse_θΨ  = fill(0.0::Float64, N_SoilSelect)
            Nse_KΨ  = fill(0.0::Float64, N_SoilSelect)
            Rmse    = fill(0.0::Float64, N_SoilSelect)
            Rmse_θΨ = fill(0.0::Float64, N_SoilSelect)
            Rmse_KΨ = fill(0.0::Float64, N_SoilSelect)

				hydroOther = HYDRO_OTHER(Nse, Nse_θΨ, Nse_KΨ, Rmse, Rmse_θΨ, Rmse_KΨ)	
				return hydroOther
			end  # function: HYDRO_OTHER
end  # module: hydroStruct

# ............................................................