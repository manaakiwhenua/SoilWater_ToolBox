# =============================================================
#		MODULE: cst
# =============================================================
module cst
   const Day_2_Second  = 60.0 * 60.0 * 24.0
   const Second_2_Day  = 1.0 / (60.0 * 60.0 * 24.0)
   const Kpa_2_Mm      = 101.97162 # conversion factor from kPa to mm H2O
   const MmS_2_CmH     = 60.0 * 60.0 / 10.0
   const MmS_2_MmH     = 60.0 * 60.0

   const Mm_2_Cm       = 0.1
   const Cm_2_Mm       = 10.0
   const Hour_2_Second = 60.0 * 60.0
   const Mm_2_kPa_Exact= 1.0 / 101.97162
   const Mm_2_kPa      = 0.01
   const kPa_2_Mm      = 100.0
   const kPa_2_Mm_Exact = 101.97162
   
   const Second_2_Hour = 1.0 / (60.0 * 60.0)
   const Y             = 14.9 #Young Laplace equation constant [mm2]
   const β             = 0.6 # Best infiltration parameter
   const KunsatModel = 119980.32407407407 # [mm s ⁻¹]

end  # module cst
# ............................................................