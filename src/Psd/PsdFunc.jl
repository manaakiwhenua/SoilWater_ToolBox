module psdFunc
	import ..option

	export PSD_MODEL
	import BlackBoxOptim
		
	# =========================================
   	#       PSD MODELS
	# ========================================
		function PSD_MODEL(iZ, Psd, ∑Psd, Rpart, N_Psd, θs_Psd, θr_Psd, paramPsd)
 		
			if option.psd.Model == :IMP
				# Correction for the small PSD
				Psd, ∑Psd = imp.SUBCLAY_CORRECTION(∑Psd, paramPsd.Subclay[iZ], N_Psd)  		

				# Computing ξ2 from Psd
					ξ2 = imp.∑PSD_2_ξ2(∑Psd[paramPsd.∑Psd_2_ξ2_Size[iZ]]; ∑Psd_2_ξ2_β1=paramPsd.∑Psd_2_ξ2_β1[iZ], ∑Psd_2_ξ2_β2=paramPsd.∑Psd_2_ξ2_β2[iZ])

					
				# Computing θ from Psd
					θ_Rpart = imp.RPART_2_θ(θs_Psd, θr_Psd, Psd, Rpart, N_Psd, paramPsd.ξ1[iZ], ξ2) 
					
				# Computing Ψ from Psd
					Ψ_Rpart = imp.RPART_2_ΨRPART(Rpart, N_Psd) 	# Computing  from Psd
	
			elseif option.psd.Model == :Chang2019Model
				# Computing θ from Psd
					θ_Rpart = chang.RPART_2_θ(θs_Psd, Psd, Rpart, N_Psd, paramPsd.ξ1[iZ]) 
				
				# Computing Ψ from Psd
					Ψ_Rpart = chang.RPART_2_ΨRPART(Rpart, N_Psd)
				 									
			end # option.psd.Chang2019Model
			
			return θ_Rpart, Ψ_Rpart
		end # function PSD_MODEL

	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>
	# <>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>	


	# =============================================================
	#		MODULE: imp
	# =============================================================
	module imp
		import ...cst, ...param, ...psdInitialize
		export ∑PSD_2_ξ2, SUBCLAY_CORRECTION, INTERGRANULARMIXING, RPART_2_ΨRPART

		# =========================================
		#      Rpart -> Ψ_Rpart
		# =========================================
			function RPART_2_ΨRPART(Rpart, N_Psd) 
				Ψ_Rpart = zeros(Float64, N_Psd)
				
				# It is to be noted that
				Rpart_Max = Rpart[N_Psd]
				Rpart_Min = Rpart[1]
				
				return Ψ_Rpart =  param.psd.imp.Ψ_Max .* ( ( (cst.Y  ./ Rpart[1:N_Psd]) .- (cst.Y ./ Rpart_Max) ) ./ ((cst.Y  ./ Rpart_Min) - (cst.Y  ./ Rpart_Max)) ) .^ param.psd.imp.λ 
			end # function RPART_2_ΨRPART


		# =========================================
		#      INTERGRANULARMIXING MODELS
		# =========================================
			function INTERGRANULARMIXING(Rpart, ξ1, ξ2)
				return IntergranularMixing = min(max(ξ1 * exp(-(Rpart ^ -ξ2)), 0.0), param.psd.imp.ξ_Max)
			end # function INTERGRANULARMIXING


		# =========================================
		#      UNIVERSAL INTERGRANULARMIXING MODEL
		# =========================================
			function ∑PSD_2_ξ2(∑Psd; ∑Psd_2_ξ2_β1=param.psd.imp.∑Psd_2_ξ2_β1, ∑Psd_2_ξ2_β2=param.psd.imp.∑Psd_2_ξ2_β2)
				return ξ2 = min(∑Psd_2_ξ2_β1 * exp(∑Psd_2_ξ2_β2 * ∑Psd), param.psd.imp.ξ2_Max)   # ξ2 = ∑Psd_2_ξ2_β1 + (∑Psd_2_ξ2_β2 * ∑Psd) 
			end # function ∑PSD_2_ξ2

			function MODEL_ξ1(P_ξ1=param.psd.imp.P_ξ1)
				return ξ1 = P_ξ1
			end # function MODEL_ξ1


		# =========================================
		#       Rpart -> θ
		# =========================================
			function RPART_2_θ(θs, θr_Psd, Psd, Rpart, N_Psd, ξ1, ξ2)
				θ_Rpart = zeros(Float64, N_Psd)

				# Computing the divisor
					∑θRpart = Psd[1] / (Rpart[1] ^ INTERGRANULARMIXING(Rpart[1], ξ1, ξ2))
					for iRpart=2:N_Psd
						∑θRpart +=  Psd[iRpart] / (Rpart[iRpart] ^ INTERGRANULARMIXING(Rpart[iRpart], ξ1, ξ2))
					end
			
				# Computing the dividend
					θ_Rpart[1] =  Psd[1] / (Rpart[1] ^ INTERGRANULARMIXING(Rpart[1], ξ1, ξ2))
					for iRpart=2:N_Psd
						θ_Rpart[iRpart] = θ_Rpart[iRpart-1] + Psd[iRpart] / (Rpart[iRpart] ^ INTERGRANULARMIXING(Rpart[iRpart], ξ1, ξ2))
					end

				# Computing θ_Rpart
					for iRpart=1:N_Psd
						θ_Rpart[iRpart] =  (θs - θr_Psd) * (θ_Rpart[iRpart] / ∑θRpart) + θr_Psd
					end

				return θ_Rpart
			end # function RPART_2_θ



		# =========================================
		#          Subclay -> ∑PSD, PSD
		# =========================================
			function SUBCLAY_CORRECTION(∑Psd, Subclay, N_Psd)
				# Correction for the small PSD
				# Subclay = 1.0 # no subclay correction applied
				∑Psd[1] = ∑Psd[1] * Subclay
				Psd = psdInitialize.∑PSD_2_PSD(∑Psd[1:N_Psd], N_Psd)
				return Psd, ∑Psd
			end # Subclay


		# =========================================
		#          PSD -> ∑PSD
		# =========================================
			function PSD_2_∑PSD(Psd, N_Psd)
				∑Psd = zeros(Float64, N_Psd)
				∑Psd[1] = Psd[1]
				 for iRpart = 2:N_Psd
					∑Psd[iRpart] = ∑Psd[iRpart-1] + Psd[iRpart]
				end
				return ∑Psd
			end # function PSD_2_∑PSD

			# =========================================
			#      Rpart -> r_pore FOR NEXT!!!!     r_pore = Rpart.*(((1.0./(θs.*ρp)).* ∑psd.*(Rpart).^(-ξ)).^(3-ξ))
			# =========================================
			# =========================================
			#      r_pore -> Ψ_Rpart FOR NEXT!!!!   Ψ_Rpart = Y ./ r_pore      #Young Laplace
			# =========================================

	end  # module: imp
	# ............................................................


	# =============================================================
	#		MODULE: Chang et al., 2019
	# =============================================================
	module chang
		import ...cst, ...param

		# ==============================================
		#      Rpart -> Ψ_Rpart  from Chang et al., 2019
		# ==============================================
			function RPART_2_ΨRPART(Rpart, N_Psd) 
				Ψ_Rpart = zeros(Float64, N_Psd)
				return Ψ_Rpart =  cst.Y ./ (0.3 .* Rpart[1:N_Psd]) 
			end # function RPART_2_ΨRPART_Chang


		# ==============================================
		#        Rpart -> θ  from Chang et al., 2019
		# ==============================================
			function RPART_2_θ(θs, Psd, Rpart, N_Psd, β)
				θ_Rpart = zeros(Float64, N_Psd)
				δ = zeros(Float64, N_Psd)
				∑ = zeros(Float64, N_Psd)

				Clay = Psd[1]
				∑[1] = Clay ^ β
				 for iRpart = 2:N_Psd
					δ[iRpart] = Psd[iRpart] / sum(Psd[2:N_Psd])
					∑[iRpart] = ∑[iRpart-1] + Psd[iRpart] - ((Clay ^ β) - Clay) * δ[iRpart]
				end

				θ_Rpart[1] = θs * ∑[1] 
				 for iRpart = 2:N_Psd
					θ_Rpart[iRpart] = θs * ∑[iRpart] 
				end
				return θ_Rpart
			end # function RPART_2_θ_Chang
	
	end  # module chang
	# ............................................................
	
end # module psdFunc