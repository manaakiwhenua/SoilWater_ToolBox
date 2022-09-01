# =============================================================
#		MODULE: residual
# =============================================================
module residual

	import ..flux, ..wrc, ..pond, ..kunsat
	import ForwardDiff: derivative
	export RESIDUAL_DIFF, RESIDUAL, ∂RESIDUAL∂Ψ, ∂R∂Ψ_FORWARDDIFF, ∂R∂Ψ▽_FORWARDDIFF, ∂R∂Ψ△_FORWARDDIFF

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	# 		FUNCTION : RESIDUAL_DIFF DERIVATIVE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RESIDUAL(discret, hydro, iT::Int64, iZ::Int64, N_iZ::Int64, Q, Residual, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ, Ψ_Max, Ψ_Min)
			if iZ==1
				Q[iT,1] = flux.Q!(discret, hydro, 1, iT, N_iZ, ΔHpond, ΔPr, ΔT, Ψ[iT,1], Ψ[iT,1])
			end

			# if iZ == N_iZ
			# 	Q[iT,iZ+1] = flux.Q!(discret, hydro, iZ+1, iT, N_iZ, ΔHpond, ΔPr, ΔT, Ψ[iT, min(iZ+1,N_iZ)], Ψ[iT,iZ])
			
			# 	θ[iT,iZ] = wrc.Ψ_2_θDual(Ψ[iT,iZ], iZ, hydro)

			# 	Qmax = ( - ΔSink[iT,iZ] - discret.ΔZ[iZ] * ((hydro.θs[iZ] - θ[iT-1,iZ]) - hydro.So[iZ] * (eps()- Ψ[iT-1,iZ]) * (θ[iT,iZ] / hydro.θs[iZ])))/ ΔT[iT ]+ Q[iT,iZ]

			# 	Q[iT,iZ+1] = min(max(Qmax, 0.0), Q[iT,iZ+1])

			# else
				Q[iT,iZ+1] = flux.Q!(discret, hydro, iZ+1, iT, N_iZ, ΔHpond, ΔPr, ΔT, Ψ[iT, min(iZ+1,N_iZ)], Ψ[iT,iZ])

				θ[iT,iZ] = wrc.Ψ_2_θDual(Ψ[iT,iZ], iZ, hydro)
		#   end
			
			Residual[iZ] = discret.ΔZ[iZ] * ((θ[iT,iZ] - θ[iT-1,iZ]) - hydro.So[iZ] * (Ψ[iT,iZ]- Ψ[iT-1,iZ]) * (θ[iT,iZ] / hydro.θs[iZ])) - ΔT[iT] * (Q[iT,iZ] - Q[iT,iZ+1]) + ΔSink[iT,iZ]

			return Q, Residual, θ
		end  # function: RESIDUAL_DIFF



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RESIDUAL_DIFF DERIVATIVE
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RESIDUAL_DIFF(Flag_NoConverge::Bool, discret, hydro, iT::Int64, iZ::Int64, N_iZ::Int64, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ▲, Ψ₀, Ψbest_, Ψbest▲, Ψbest▼, Ψ_, Ψ▼, Ψ_Max)

			if !Flag_NoConverge
				# Q[iT,iZ] format for ForwardDiff
					Q₁ = flux.Q!(discret, hydro, iZ, iT, N_iZ, ΔHpond, ΔPr, ΔT, Ψ_, Ψ▲)

				# Q[iT,iZ+1] format for ForwardDiff
					Q₂ = flux.Q!(discret, hydro, iZ+1, iT, N_iZ, ΔHpond, ΔPr, ΔT, Ψ▼, Ψ_)		
			else 
				# Q[iT,iZ] format for ForwardDiff
					Q₁ = flux.Q!(discret, hydro, iZ, iT, N_iZ, ΔHpond, ΔPr, ΔT, Ψbest_, Ψbest▲)

				# Q[iT,iZ+1] format for ForwardDiff
					Q₂ = flux.Q!(discret, hydro, iZ+1, iT, N_iZ, ΔHpond, ΔPr, ΔT, Ψbest▼, Ψbest_)
			end

			# θ[iT,iZ] format for ForwardDiff
				θ₂ = wrc.kg.Ψ_2_θDual(Ψ_, iZ, hydro)

		return Residual₂ = discret.ΔZ[iZ] * ((θ₂ - θ[iT-1,iZ]) - hydro.So[iZ] * (Ψ_ - Ψ₀) * (θ[iT,iZ] / hydro.θs[iZ])) - ΔT[iT] * (Q₁ - Q₂) + ΔSink[iT,iZ] 
		end  # function: RESIDUAL_DIFF



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ∂R∂Ψ_NUMERICAL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ∂R∂Ψ_FORWARDDIFF(Flag_NoConverge::Bool, discret, hydro, iT::Int64, iZ::Int64, N_iZ::Int64, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ▲, Ψ₀, Ψbest_, Ψbest▲, Ψbest▼, Ψ_, Ψ▼, Ψ_Max)	

			ψ = Ψ_

			∂R∂Ψ_Func(ψ) = RESIDUAL_DIFF(Flag_NoConverge, discret, hydro, iT, iZ, N_iZ, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ▲, Ψ₀, Ψbest_, Ψbest▲, Ψbest▼, max(ψ,0.0), Ψ▼, Ψ_Max)[1]
			
			∂R∂Ψ_Derivative_1 = ψ -> derivative(∂R∂Ψ_Func, ψ)	

			return ∂R∂Ψ = ∂R∂Ψ_Derivative_1(ψ)
		end # function: ∂RESIDUAL∂Ψ_NUMERICAL


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ∂R∂Ψ▽_NUMERICAL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ∂R∂Ψ▽_FORWARDDIFF(Flag_NoConverge::Bool, discret, hydro, iT::Int64, iZ::Int64, N_iZ::Int64, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ▲, Ψ₀, Ψbest_, Ψbest▲, Ψbest▼, Ψ_, Ψ▼, Ψ_Max)

			if iZ ≤ N_iZ-1
				ψ▼ = Ψ▼

				∂R∂Ψ_Func(ψ▼) = RESIDUAL_DIFF(Flag_NoConverge, discret, hydro, iT, iZ, N_iZ, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ▲, Ψ₀, Ψbest_, Ψbest▲, Ψbest▼, Ψ_, max(ψ▼,0.0), Ψ_Max)[1]

				∂R∂Ψ_Derivative_1 = ψ▼ -> derivative(∂R∂Ψ_Func, ψ▼)			
				
				return ∂R∂Ψ▽ = ∂R∂Ψ_Derivative_1(ψ▼)
			else
				return ∂R∂Ψ▽ = 0.0
			end
		end # function: ∂RESIDUAL∂Ψ_NUMERICAL


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ∂R∂Ψ▽_NUMERICAL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ∂R∂Ψ△_FORWARDDIFF(Flag_NoConverge::Bool, discret, hydro, iT::Int64, iZ::Int64, N_iZ::Int64, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ▲, Ψ₀, Ψbest_, Ψbest▲, Ψbest▼, Ψ_, Ψ▼, Ψ_Max)

			if iZ ≥ 2
				ψ▲ = Ψ▲

				∂R∂Ψ_Func(ψ▲) =  RESIDUAL_DIFF(Flag_NoConverge, discret, hydro, iT, iZ, N_iZ, ΔHpond, ΔPr, ΔSink, ΔT, θ, max(ψ▲,0.0), Ψ₀, Ψbest_,Ψbest▲, Ψbest▼, Ψ_, Ψ▼, Ψ_Max)[1]
				
				∂R∂Ψ_Derivative_1 = ψ▲ -> derivative(∂R∂Ψ_Func, ψ▲)			
				
				return ∂R∂Ψ△ = ∂R∂Ψ_Derivative_1(ψ▲)
			else
				return ∂R∂Ψ△ = 0.0
			end
		end # function: ∂RESIDUAL∂Ψ_NUMERICAL


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION2 : ∂∂R∂Ψ_NUMERICAL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ∂∂R∂Ψ_FORWARDDIFF(Flag_NoConverge::Bool, discret, hydro, iT::Int64, iZ::Int64, N_iZ::Int64, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ▲, Ψ₀, Ψbest_, Ψbest▲, Ψbest▼, Ψ_, Ψ▼, Ψ_Max)	
			ψ = Ψ_

			∂R∂Ψ_Func(ψ) = RESIDUAL_DIFF(Flag_NoConverge, discret, hydro, iT, iZ, N_iZ, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ▲, Ψ₀, Ψbest_, Ψbest▲, Ψbest▼, max(ψ,0.0), Ψ▼, Ψ_Max)[1]
			
			∂R∂Ψ_Derivative_1 = ψ -> derivative(∂R∂Ψ_Func, ψ)	

			∂R∂Ψ_Derivative_2 = ψ -> derivative(∂R∂Ψ_Derivative_1, ψ)	
			
			∂R∂Ψ = ∂R∂Ψ_Derivative_1(ψ)

			∂∂R∂Ψ = ∂R∂Ψ_Derivative_2(ψ)

			return ∂R∂Ψ, ∂∂R∂Ψ
		end # function: ∂RESIDUAL∂Ψ_NUMERICAL


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION2 : ∂R∂Ψ▽_NUMERICAL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ∂∂R∂Ψ▽_FORWARDDIFF(Flag_NoConverge::Bool, discret, hydro, iT::Int64, iZ::Int64, N_iZ::Int64, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ▲, Ψ₀, Ψbest_, Ψbest▲, Ψbest▼, Ψ_, Ψ▼, Ψ_Max)
			ψ▼ = Ψ▼

			∂R∂Ψ_Func(ψ▼) = RESIDUAL_DIFF(Flag_NoConverge, discret, hydro, iT, iZ, N_iZ, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ▲, Ψ₀, Ψbest_, Ψbest▲, Ψbest▼, Ψ_, max(ψ▼,0.0), Ψ_Max)[1]
			
			∂R∂Ψ_Derivative_1 = ψ▼ -> derivative(∂R∂Ψ_Func, ψ▼)
			
			∂R∂Ψ_Derivative_2 = ψ▼ -> derivative(∂R∂Ψ_Derivative_1 , ψ▼)	

			return ∂∂R∂Ψ▽2 = ∂R∂Ψ_Derivative_2(ψ▼)
		end # function: ∂RESIDUAL∂Ψ_NUMERICAL


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION2 : ∂R∂Ψ▽_NUMERICAL
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ∂∂R∂Ψ△_FORWARDDIFF(Flag_NoConverge::Bool, discret, hydro, iT::Int64, iZ::Int64, N_iZ::Int64, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ▲, Ψ₀, Ψbest_, Ψbest▲, Ψbest▼, Ψ_, Ψ▼, Ψ_Max)
			ψ▲ = Ψ▲

			∂R∂Ψ_Func(ψ▲) = RESIDUAL_DIFF(Flag_NoConverge, discret, hydro, iT, iZ, N_iZ, ΔHpond, ΔPr, ΔSink, ΔT, θ, max(ψ▲,0.0), Ψ₀, Ψbest_,Ψbest▲, Ψbest▼, Ψ_, Ψ▼, Ψ_Max)[1]
			
			∂R∂Ψ_Derivative_1 = ψ▲ -> derivative(∂R∂Ψ_Func, ψ▲)
			
			∂R∂Ψ_Derivative_2 = ψ▲ -> derivative(∂R∂Ψ_Derivative_1, ψ▲)		
			
			# ∂R∂Ψ△1 = ∂R∂Ψ_Derivative_2(ψ▲)

			return ∂∂R∂Ψ△2 = ∂R∂Ψ_Derivative_2(ψ▲)
		end # function: ∂RESIDUAL∂Ψ_NUMERICAL


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ∂RESIDUAL∂Ψ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function ∂RESIDUAL∂Ψ(∂K∂Ψ, discret, hydro, iT::Int64, iZ::Int64,  N_iZ::Int64, ΔT, θ, Ψ)

			# K_Aver[iZ+1] = flux.K_AVER!(discret, hydro, iZ+1, N_iZ, Ψ[iT, min(iZ+1,N_iZ)],  Ψ[iT, iZ])

			Sw = hydro.So[iZ] / hydro.θs[iZ]

			∂K∂Ψ[iZ] = kunsat.∂K∂Ψ(Ψ[iT,iZ], iZ, hydro)

			∂θ∂Ψ = wrc.∂θ∂Ψ(Ψ[iT,iZ], iZ, hydro)

			∂Q∂Ψ = flux.∂q∂Ψ.∂Q∂Ψ(∂K∂Ψ, discret, hydro, iT, iZ,  N_iZ, Ψ)

			∂Q∂Ψ△ = flux.∂q∂Ψ.∂Q∂Ψ△(∂K∂Ψ, discret, hydro, iT, iZ,  N_iZ, Ψ)

			∂Q▽∂Ψ = flux.∂q∂Ψ.∂Q▽∂Ψ(∂K∂Ψ, discret, hydro, iT, iZ+1,  N_iZ, Ψ)

			∂Q▽∂Ψ▽ = flux.∂q∂Ψ.∂Q▽∂Ψ▽(∂K∂Ψ, discret, hydro, iT, iZ,  N_iZ, Ψ)
		
			if iZ ≥ 2
				∂R∂Ψ△ = - ΔT[iT] * ∂Q∂Ψ△
			else
				∂R∂Ψ△ = 0.0
			end

			∂R∂Ψ = discret.ΔZ[iZ] * (∂θ∂Ψ * (1.0 - Sw * (Ψ[iT,iZ] - Ψ[iT-1,iZ]) ) - Sw * θ[iT,iZ]) - ΔT[iT] * (∂Q∂Ψ - ∂Q▽∂Ψ)
		
			if iZ ≤ N_iZ-1
				∂R∂Ψ▽ = ΔT[iT] * ∂Q▽∂Ψ▽
			else
				∂R∂Ψ▽ = 0.0
			end

			return ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽
		end  
	
end  # module: residual
# ............................................................