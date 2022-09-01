module interpolate
	export ∑_2_Δ, POINTS_2_SlopeIntercept, INTERPOLATE_2D_LOOP

	"""
	∑RAIN_2_ΔPr(∑Pr, ∑Pr_Climate, ∑T_Climate, ∑T, iT, iT_Pr, N_Climate; FlagForwardTime=true)
	This is used for other variables, we give a example for Pr
	"""
	function ∑_2_Δ(∑X_Past, ∑X_Climate, ∑T, ∑T_Climate, iT_X, N_Climate, Flag_ReRun, iT)

		# Moving backwards if we need to rerun
			if Flag_ReRun && iT_X ≥ 3
				iT_X -= 1
			end

		# Determening if we should increase iT_X
			FlagBreak = false
			while !(FlagBreak)
				if (∑T_Climate[iT_X-1]≤ ∑T[iT] ≤ ∑T_Climate[iT_X]) || (iT_X == N_Climate) 
					FlagBreak = true
					break
				else 
					iT_X += 1
					FlagBreak = false
				end # if
			end # while

		# Building a regression line which passes from POINT1(∑T_Climate[iT_X], ∑Pr_Climate[iT_Pr]) and POINT2: (∑T_Climate[iT_Pr+1], ∑Pr_Climate[iT_Pr+1])
			Slope, Intercept = POINTS_2_SlopeIntercept(∑T_Climate[iT_X-1], ∑X_Climate[iT_X-1], ∑T_Climate[iT_X], ∑X_Climate[iT_X])

			∑X = Slope * ∑T[iT] + Intercept

		# Xecipitation [mm /  ΔTconst]
			ΔX = ∑X - ∑X_Past
		
		return ∑X, ΔX, iT_X

	end # function ∑RAIN_2_ΔX


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : POINTS_2_SlopeIntercept
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	"""POINTS_2_SlopeIntercept
	From Point1 [X1, Y1] and point2 [X2, Y2] compute Y = Slope.X + Intercept
	"""
		function POINTS_2_SlopeIntercept(X1, Y1, X2, Y2)
			Slope = (Y2 - Y1) / (X2 - X1)
			Intercept = (Y1 * X2 - X1 * Y2) / (X2 - X1)
			return Slope, Intercept
		end # POINTS_2_SlopeIntercept


	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INTERPOLATE_2D_LOOP
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INTERPOLATE_2D_LOOP(∑T, ∑T_Reduced, NiT_Reduced, N_iT, N_iZ, X_Reduced, X)
			iT_X = 2
			iCount_TooEarly = 0 

			for iT=1:NiT_Reduced
				FlagBreak = false
				Flag_TooEarly = false
				while !(FlagBreak)
					if (∑T[iT_X-1] ≤ ∑T_Reduced[iT] ≤ ∑T[iT_X]) || (iT_X == N_iT) 
						FlagBreak = true
						break
					elseif ∑T_Reduced[iT] < ∑T[iT_X-1]
						Flag_TooEarly = true
						iCount_TooEarly += 1
						break			
					else 
						iT_X += 1
						FlagBreak = false
					end # if
				end # while

				# Building a regression line which passes from POINT1(∑T_Climate[iT_X], ∑Pr_Climate[iT_Pr]) and POINT2: (∑T_Climate[iT_Pr+1], ∑Pr_Climate[iT_Pr+1])
				if !Flag_TooEarly
					for iZ = 1:N_iZ
						Slope, Intercept = interpolate.POINTS_2_SlopeIntercept(∑T[iT_X-1], X[iT_X-1,iZ], ∑T[iT_X], X[iT_X,iZ])

						X_Reduced[iT,iZ] = Slope * ∑T_Reduced[iT] + Intercept
					end # for iZ = 1:N_iZ
				else
					for iZ = 1:N_iZ
						X_Reduced[iT,iZ] = X[iT_X,iZ] # TODO problem of shifting of 1 day
					end
				end
			end # for: iT=1:obsθ.N_iT

			# TODO to be checked
			if iCount_TooEarly ≥ 1
				for iZ = 1:N_iZ
					X_Reduced[1:NiT_Reduced-1, iZ] = X_Reduced[2:NiT_Reduced, iZ]
					X_Reduced[NiT_Reduced, iZ] = X[NiT_Reduced,iZ]
				end

			end
		
			
		return X_Reduced
	end  # function: θINTERPOLATION

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : INTERPOLATE_2D_LOOP
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		function INTERPOLATE_1D_LOOP(∑T, ∑T_Reduced, NiT_Reduced, N_iT, X_Reduced, X)
			iT_X = 2
			for iT=1:NiT_Reduced
		
				FlagBreak = false
				while !(FlagBreak)
					if (∑T[iT_X-1] - eps(10.0) ≤ ∑T_Reduced[iT] ≤ ∑T[iT_X] + eps(10.0)) || (iT_X == N_iT) 
						FlagBreak = true
						break
					else 
						iT_X += 1
						FlagBreak = false
					end # if
				end # while

				# Building a regression line which passes from POINT1(∑T_Climate[iT_X], ∑Pr_Climate[iT_Pr]) and POINT2: (∑T_Climate[iT_Pr+1], ∑Pr_Climate[iT_Pr+1])
				Slope, Intercept = interpolate.POINTS_2_SlopeIntercept(∑T[iT_X-1], X[iT_X-1], ∑T[iT_X], X[iT_X])

				X_Reduced[iT] = Slope * ∑T_Reduced[iT] + Intercept
			
			end # for: iT=1:obsθ.N_iT
		
		return X_Reduced
	end  # function: θINTERPOLATION

end # module interpolate