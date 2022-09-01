# =============================================================
#		module: θINITIAL()
# =============================================================
module θini
   import ..wrc
   export θINI_TOP

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : θINI_TOP
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   """ 
      θINI_TOP(hydro, N_iZ, θ)
   Compute θ[t=1:1:N_iZ] when only θ[1,1] is known
   It assumes that Se[t=1:1:N_iZ] is constant through the soil profile"""
      function θINI_TOP(hydro, N_iZ, θ)
         Se = wrc.θ_2_Se(θ[1,1], 1, hydro)

         for iZ=2:N_iZ
            θ[1,iZ] = wrc.Se_2_θ(Se, 1, hydro) 
         end
      return θ
      end  # function: θINI_TOP

end  # module: θini
# ............................................................