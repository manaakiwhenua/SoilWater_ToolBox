module timeTransSteady
	import ..stats, ..param
	export ∑INFIlT_2_TIMETRANSSTEADY

	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#		FUNCTION : ∑INFIlT_2_TIMETRANSSTEADY
	#		Determening from the data when the transition between transit and steady occures
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~		
	function  ∑INFIlT_2_TIMETRANSSTEADY(T, N_SoilSelect, N_Infilt, infiltOutput, ∑Infilt_Obs; N_LastInfiltPoint=4, CorrectLinear=true) 

		# FOR EVERY SOIL
		for iZ=1:N_SoilSelect

			Flag_Break = false

			# Want at least 3 remaining points
				if N_Infilt[iZ] - N_LastInfiltPoint <= 3
					# Must at least leave 3 points
					N_LastInfiltPoint = max(N_Infilt[iZ] - 3, 2)	
				end	

			# Determine when it is no longer linear
				iStart = N_Infilt[iZ] - N_LastInfiltPoint + 1
				iEnd = N_Infilt[iZ]

				Intercept, Slope = stats.LINEAR_REGRESSION(T[iZ,iStart:iEnd], ∑Infilt_Obs[iZ,iStart:iEnd])

			# Starting from the last soils
				for i =1:N_Infilt[iZ] - N_LastInfiltPoint - 1
					iModel = N_Infilt[iZ] - N_LastInfiltPoint - i
					# iEnd = N_Infilt[iZ] - N_LastInfiltPoint

					# Determine the linear regression
						∑Infilt_Model = T[iZ,iModel] * Slope + Intercept

					#Determine if enough points for the linear regression since it must monotically decrease
						if ∑Infilt_Model > ∑Infilt_Obs[iZ,iModel] && CorrectLinear == true
								# Recompute the slope and intercept
							Intercept, Slope = stats.LINEAR_REGRESSION(T[iZ,iModel:iEnd], ∑Infilt_Obs[iZ,iModel:iEnd])
							∑Infilt_Model = T[iZ,iModel] * Slope + Intercept
						end

					# Compute the error of slope of not fitting the linear steady equation
						ΔSlope_Err = abs(∑Infilt_Model - ∑Infilt_Obs[iZ,iModel]) / (T[iZ,iModel+1]-T[iZ,iModel])
						ΔSlope_Err = rad2deg(atan(abs(ΔSlope_Err)))
					
					if (ΔSlope_Err >= param.infilt.ΔSlope_Err_SteadyState || iModel<=3) && Flag_Break == false

						Flag_Break = true

						iModel = max(iModel - 1, 3)

						infiltOutput.iT_TransSteady_Data[iZ] = iModel
						
						infiltOutput.T_TransSteady_Data[iZ] = T[iZ,iModel]

						break # To speed up
					end # 	if i-1 >= 1
					
				end # for i = N_Infilt[iZ] - N_LastInfiltPoint:-1:1

		end # for iZ=1

		return infiltOutput

	end # function: INFIlTOBS_2_iTIME_TRANS_STEADy

end  # macro timeTransSteady

