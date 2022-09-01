# =============================================================
#		MODULE: residual
# =============================================================
module richard
	import ..timeStep, ..flux, ..kunsat, ..option, ..param, ..pond, ..residual, ..wrc
	using LinearAlgebra

	export RICHARD_ITERATION

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RICHARD_ITERATION
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RICHARD_ITERATION(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, Count_ReRun, discret, Flag_NoConverge::Bool, hydro, iNonConverge::Int64, iT::Int64, Iter_CountTotal::Int64, N_iZ::Int64, Q, Residual, Sorptivity::Float64, ΔHpond, ΔΨmax, ΔPr, ΔSink, ΔT, Δθ_Max::Float64, θ, Ψ, Ψ_Max, Ψ_Min, Ψbest)

			# INITIALIZING
				for iZ=1:N_iZ
					Ψ[iT,iZ] = Ψbest[iZ]
				end
		
				# Residual_Max_Best it should be improved compared to Ψ[iT,iZ] = Ψ[iT-1,iZ]
				~, ~, ~, ~, Residual, ~, ~ = richard.RICHARD(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, discret, Flag_NoConverge, hydro, iT, N_iZ, Q, Residual, Sorptivity, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ, Ψ_Max, Ψ_Min, Ψbest)

				# Averaging the Residuals, depending on method
					Residual_Max_Best = RESIDUAL_MAX(discret, iT, N_iZ, Residual, ΔT)
	
			# ITTERATION
			iTer = 0::Int64
			while iTer ≤ param.hyPix.N_Iter	
				iTer += 1
				Iter_CountTotal += 1 # Counting the iterations
		
				∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, Q, Residual, ΔHpond, θ = richard.RICHARD(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, discret, Flag_NoConverge, hydro, iT, N_iZ, Q, Residual, Sorptivity, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ, Ψ_Max, Ψ_Min, Ψbest)

				Ψ = SOLVING_TRIAGONAL_MATRIX(∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, iT, iTer, N_iZ, Residual, ΔΨmax, Ψ, Ψ_Max, Ψ_Min)

				# Averaging the Residuals, depending on method
					Residual_Max = RESIDUAL_MAX(discret, iT, N_iZ, Residual, ΔT)

				# Determine if itteration made improvement
					if Residual_Max < Residual_Max_Best	
						for iZ=1:N_iZ
							Ψbest[iZ] = Ψ[iT,iZ]
						end
						Residual_Max_Best = Residual_Max
					end # Residual_Max < Residual_Max_Best 	

				# Did we achieve the goals
				if Residual_Max < param.hyPix.WaterBalanceResidual_Max
					break # Move out the loop
				end  # if: Residual
			end # while: iTer ======================

			# Making sure we get the best if convergence fails
			if iTer ≥ param.hyPix.N_Iter + 1
				iNonConverge += 1

				# if non converge compute Q(Ψbest)
				# if option.hyPix.NoConverge_Ψbest
				# 	Flag_NoConverge = true
				# end

				# Put the best values
				for iZ=1:N_iZ
					Ψ[iT,iZ] = Ψbest[iZ]
				end
			else
				Flag_NoConverge = false
			end #  iTer == param.hyPix.N_Iter

			if option.hyPix.Qbottom_Correction
				Q[iT,N_iZ+1] = max(( - ΔSink[iT,N_iZ] - discret.ΔZ[N_iZ] * ((θ[iT,N_iZ] - θ[iT-1,N_iZ]) - hydro.So[N_iZ] * (Ψ[iT,N_iZ]- Ψ[iT-1,N_iZ]) * (θ[iT,N_iZ] / hydro.θs[N_iZ]))) / ΔT[iT] + Q[iT,N_iZ], 0.0)
			end

			# Determine if the simulation is going to rerun with a different time step
				Count_ReRun, Flag_ReRun, ΔT = RERUN_HYPIX(Count_ReRun, discret, Flag_NoConverge, hydro, iT, N_iZ, Q, ΔHpond, ΔΨmax, ΔPr, ΔSink, ΔT, θ, Ψ)

			return Count_ReRun, Flag_NoConverge, Flag_ReRun, iNonConverge, Iter_CountTotal, Q, ΔHpond, ΔT, θ, Ψ
		end  # function: RICHARD_SOLVING



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RICHARD
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RICHARD(∂K∂Ψ, ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, discret, Flag_NoConverge::Bool, hydro, iT::Int64, N_iZ::Int64, Q, Residual, Sorptivity, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ, Ψ_Max, Ψ_Min, Ψbest)

			ΔHpond = pond.PONDING_SORPTIVITY!(discret, hydro, iT, Sorptivity, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ)

			# ∂R∂Ψ2 = fill(0.0, N_iZ)
			# ∂R∂Ψ▽2 = fill(0.0, N_iZ)
			# ∂R∂Ψ△2 =  fill(0.0, N_iZ)

			for iZ=1:N_iZ		
				Q, Residual, θ = residual.RESIDUAL(discret, hydro, iT, iZ, N_iZ, Q, Residual, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ, Ψ_Max, Ψ_Min)

				if option.hyPix.∂R∂Ψ_Numerical				
					∂R∂Ψ[iZ] = residual.∂R∂Ψ_FORWARDDIFF(Flag_NoConverge, discret, hydro, iT, iZ, N_iZ, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ[iT, max(iZ-1,1)], Ψ[iT-1, iZ], Ψ[iT-1,iZ], Ψ[iT-1,max(iZ-1,1)], Ψ[iT-1, min(iZ+1,N_iZ)], Ψ[iT, iZ], Ψ[iT, min(iZ+1,N_iZ)], Ψ_Max)

					∂R∂Ψ▽[iZ]  = residual.∂R∂Ψ▽_FORWARDDIFF(Flag_NoConverge, discret, hydro, iT, iZ, N_iZ, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ[iT, max(iZ-1,1)], Ψ[iT-1, iZ], Ψ[iT-1,iZ], Ψ[iT-1,max(iZ-1,1)], Ψ[iT-1, min(iZ+1,N_iZ)], Ψ[iT, iZ], Ψ[iT, min(iZ+1,N_iZ)], Ψ_Max)
		
					∂R∂Ψ△[iZ]  = residual.∂R∂Ψ△_FORWARDDIFF(Flag_NoConverge, discret, hydro, iT, iZ, N_iZ, ΔHpond, ΔPr, ΔSink, ΔT, θ, Ψ[iT, max(iZ-1,1)], Ψ[iT-1, iZ], Ψ[iT-1,iZ], Ψ[iT-1,max(iZ-1,1)], Ψ[iT-1, min(iZ+1,N_iZ)], Ψ[iT, iZ], Ψ[iT, min(iZ+1,N_iZ)], Ψ_Max)
				else		
					∂R∂Ψ[iZ], ∂R∂Ψ△[iZ], ∂R∂Ψ▽[iZ] = residual.∂RESIDUAL∂Ψ(∂K∂Ψ, discret, hydro, iT, iZ, N_iZ, ΔT, θ, Ψ)

				end # if option.hyPix.∂R∂Ψ_Numerical

			end #for iZ= 1:N_iZ

			# println("∂R∂Ψ=" , ∂R∂Ψ[N_iZ]," , ", ∂R∂Ψ2[N_iZ],"\n")
			# println("∂R∂Ψ▽=", ∂R∂Ψ▽[N_iZ] .- ∂R∂Ψ▽2[N_iZ], "\n")
			# println("\n")

			return ∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, Q, Residual, ΔHpond, θ
		end # function RICHARD

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RESIDUAL_MAX
	#     Averaging the Residuals, depending on method
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RESIDUAL_MAX(discret, iT::Int64, N_iZ::Int64, Residual, ΔT)

				Residual_Norm = 0.0
				Residual_Max  = 0.0

				# Does not take into consideration the last cell which has a perfect mass balance
				for iZ=1:N_iZ
					if option.hyPix.NormMin == :Norm
						Residual_Norm += (Residual[iZ] / (ΔT[iT] * discret.ΔZ[iZ])) ^ 2
					else
						Residual_Max = max( Residual_Max, abs(Residual[iZ]) / (ΔT[iT] * discret.ΔZ[iZ]) ) 
					end
				end # for: iZ=N_iZ

				if option.hyPix.NormMin == :Norm
					Residual_Max = √(Residual_Norm / N_iZ)
				end
			
		return Residual_Max
		end  # function: RESIDUAL_MAX



	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : SOLVING_TRIAGONAL_MATRIX
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function SOLVING_TRIAGONAL_MATRIX(∂R∂Ψ, ∂R∂Ψ△, ∂R∂Ψ▽, iT::Int64, iTer::Int64, N_iZ::Int64, Residual, ΔΨmax, Ψ, Ψ_Max, Ψ_Min)
	
			Matrix_Trid = Tridiagonal(∂R∂Ψ△[2:N_iZ], ∂R∂Ψ[1:N_iZ], ∂R∂Ψ▽[1:N_iZ-1])

			# Transforming from row to column
			Residual = reshape(Residual, N_iZ, 1)

			NewtonStep = Matrix_Trid \ -Residual

			for iZ=1:N_iZ
				# Iteration k-1
					Ψ₀ = Ψ[iT,iZ]
				
				# Updating Ψ
					if isnan(NewtonStep[iZ])

						Ψ[iT,iZ] = Ψ₀
					else
						Ψ[iT,iZ] += NewtonStep[iZ]

						# Making sure it is within the feasible band 
							Ψ[iT,iZ] = min(max(Ψ[iT,iZ], eps()), 10.0^7)
							# Ψ_Max[iZ]

							Ψ[iT,iZ] = param.hyPix.NewtonStepWeaken * Ψ[iT,iZ] + (1.0 - param.hyPix.NewtonStepWeaken) * Ψ₀
					end

				# # Under ralation of Ψ Averaging to help with convergence
				# 	Ψ[iT,iZ] = param.hyPix.NewtonStepWeaken * Ψ[iT,iZ] + (1.0 - param.hyPix.NewtonStepWeaken) * Ψ₀

				# if isnan(Ψ[iT,iZ])		
				# 	println(iT," , " ,Ψ[iT,:])
				# 	println("")
				# 	@show NewtonStep
				# 	println("")
				# 	@show Residual
				# 	println(" ")
				# 	@show ∂R∂Ψ
				# 	println(" ")
				# 	@show ∂R∂Ψ△
				# 	println(" ")
				# 	@show ∂R∂Ψ▽
				# 	println(" ")
				# 	error("NaN = true")
				# end
			end # for iZ=1:N_iZ
			
		return Ψ
		end  # function: SOLVING_TRIAGONAL_MATRIX


		
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : Ψ_Constrain_K1
	# 		Making sure that the steps og NR are not too big and within the limits of Δθ_Max
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function Ψ_Constrain_K1(iT, iZ, ΔΨmax, Ψ)
			if Ψ[iT,iZ] ≤ Ψ[iT-1,iZ]		
				Ψ▽ = exp(log(Ψ[iT-1,iZ]) - ΔΨmax[iZ])

				Ψ[iT,iZ] = max(Ψ[iT,iZ], Ψ▽)

			else
				Ψ△ = exp(log(Ψ[iT-1,iZ]) + ΔΨmax[iZ])
				
				Ψ[iT,iZ] = min(Ψ[iT,iZ], Ψ△)
			end

			return Ψ
		end  # function: Ψ_Constrain_K1
		

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : RERUN_HYPIX
	# 		WITH UPDATED Ψ
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function RERUN_HYPIX(Count_ReRun::Int64, discret, Flag_NoConverge::Bool, hydro, iT::Int64, N_iZ::Int64, Q, ΔHpond, ΔΨmax, ΔPr, ΔSink, ΔT, θ, Ψ)
			# Rerun if updated ΔT is smaller compared to previously Ψ

			if option.hyPix.Flag_ReRun	&& Count_ReRun ≤ 3	

				for iZ= 1:N_iZ
					Q[iT,iZ] = flux.Q!(discret, hydro, iZ, iT, N_iZ, ΔHpond, ΔPr, ΔT, Ψ[iT,iZ], Ψ[iT,max(iZ-1,1)])
				end # for: iZ= 1:N_iZ+1

				ΔT_New, ~ = timeStep.ADAPTIVE_TIMESTEP(discret, hydro, iT, N_iZ, Q, ΔΨmax, ΔSink, θ, Ψ)

				# Rerun if the new time step is smaller that the older time step
				# if Flag_NoConverge # <>=<>=<>=<>=<>
				# 	ΔT[iT] = ΔT_New
				# 	Flag_ReRun = true
				# 	Count_ReRun += 1

				if ΔT[iT] / ΔT_New > param.hyPix.ΔT_Rerun # <>=<>=<>=<>=<>
					ΔT[iT] = ΔT_New
					Flag_ReRun = true
					Count_ReRun += 1
				
				else # <>=<>=<>=<>=<>
					Flag_ReRun = false
					Flag_NoConverge = false
					Count_ReRun = 0
				end
			else
				Flag_ReRun = false
				Flag_NoConverge = false
				Count_ReRun = 0
			end  # if: param.hyPix.ΔT_Rerun

			return Count_ReRun, Flag_ReRun, ΔT
		end  # function: RERUN_HYPIX

end # module: richard	