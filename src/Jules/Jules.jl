# =============================================================
#		module: jules joseph2
# =============================================================
module jules
   import ..option, ..param, ..path, ..tool, ..θini, ..hydroStruct, ..reading, ..tool, ..wrc, ..pathHypix
   import DelimitedFiles, Dates, CSV, Tables, NCDatasets, NetCDF
   export READ_JULES_SITES

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : START_JULES
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   function START_JULES()

      Options_θjules = "BrookCorey" # <"Texture">,<"vanGenuchten">,<"BrookCorey">

      pathHyPix = pathHypix.PATHHYPIX(1)

      println("    ~  $(pathHyPix.JulesMetadata) ~")

      Path_Climate = "D:\\DATAraw\\JULESdata\\Climate\\VCSN_Obs\\"

      Path_θ       = "D:\\DATAraw\\JULESdata\\SoilMoisture\\SoilMoisture_Site\\"

      Path_θJules  = "D:\\DATAraw\\JULESdata\\SoilMoisture_Jules\\SoilMoistureJules_Site\\"

      # Read data
         Data = DelimitedFiles.readdlm(pathHyPix.JulesMetadata, ',')
      # Read header
         Header = Data[1,1:end]
      # Remove first READ_ROW_SELECT
         Data = Data[2:end,begin:end]
      # Reading
         SiteName, NsiteName   = tool.readWrite.READ_HEADER_FAST(Data, Header, "SiteName")
         VCSNgridnumber, ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, "VCSNgridnumber")
         SiteNumber, ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, "SiteNumber")
         SoilName, ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, "SoilName-Smap")

      # Dictionary
      # SiteName_2_VCSNgridnumber::Dict{String, Int64} 
         iSite = 1
         SiteName_2_VCSNgridnumber = Dict("a"=>9999) # Initializing
         SiteName_2_SiteNumber  = Dict("a"=>9999)
         SoilName_2_SiteName = Dict("a"=>"b")
         SiteName_2_θini = Dict("a"=>9999.0) # Initializing
             
         for iSiteName in SiteName

            # Making a new path if not exist
               Path_Output =  path.Home * "//INPUT//DataHyPix//JULES//" * iSiteName
               mkpath(Path_Output) 
               
            # dictionary which correspond SiteName to VCSNgridnumber
               SiteName_2_VCSNgridnumber[iSiteName] = VCSNgridnumber[iSite]
               SiteName_2_SiteNumber[iSiteName] = SiteNumber[iSite]
               SoilName_2_SiteName[SoilName[iSite]] = iSiteName
          
            # Reading climate
               Path_Climate_Input = Path_Climate * "VCSN_obsSM_" * string(SiteName_2_VCSNgridnumber[iSiteName]) * ".csv"

               Path_Climate_Output = Path_Output * "//" * iSiteName * "_Daily_Climate_2.csv"
               
               READ_WRIITE_CLIMATE(Path_Climate_Input, Path_Climate_Output)
         
            # Reading obs θ
               Path_θ_Input = Path_θ * "sm_obs_" * string(SiteName_2_SiteNumber[iSiteName]) * ".nc"

               Path_θ_Output = Path_Output * "//" * iSiteName * "_Soilmoisture.csv"

               Path_Date_Output = Path_Output * "//" * iSiteName * "_Dates.csv"

               θᵢₙᵢ= READ_WRITE_θobs(iSiteName, Path_Date_Output, Path_θ_Input, Path_θ_Output)

               SiteName_2_θini[iSiteName] = θᵢₙᵢ
      

               # Reading Jules simulated θ
                  # Path_θjules_Input = Path_θJules * "Sta_" * string(SiteName_2_SiteNumber[iSiteName]) * "/"
   
                  # Path_θjules_Output =  Path_Output * "//" * iSiteName * "_Soilmoisture_Jules.csv"
   
                  # READ_WRITE_θJULES(Path_θjules_Input, Path_θjules_Output, Options_θjules)
            iSite += 1
         end #  for iSiteName

   return SoilName_2_SiteName,  SiteName_2_θini
   end  # function: START_JULES

   
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : READ_CLIMATE
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function READ_WRIITE_CLIMATE(Path_Input, Path_Output)

         # Read data
            Data = DelimitedFiles.readdlm(Path_Input, ',')
         # Read header
            Header = Data[1,1:end]
         # Remove first READ_ROW_SELECT
            Data = Data[2:end,begin:end]
         # Reading
            Date, N = tool.readWrite.READ_HEADER_FAST(Data, Header, "OBS_DATE")
            Rain, ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, "rain_fill")
            Pet, ~ = tool.readWrite.READ_HEADER_FAST(Data, Header, "pet_fao56")

            Year  = Dates.year.(Dates.DateTime.(Date))
            Month = Dates.month.(Dates.DateTime.(Date))
            Day   = Dates.day.(Dates.DateTime.(Date))
            Hour  = fill(9::Int64, N)
            Minute = fill(0::Int64, N)
            Second = fill(0::Int64, N)
            
            Header = ["Year";"Month";"Day";"Hour";"Minute";"Second";"PET(mm)";"Rain(mm)"]

            Output = Tables.table( [Year Month Day Hour Minute Second Pet Rain])

            CSV.write(Path_Output, Output, header=Header)	
            
         return nothing
         end  # function: READ_CLIMATE


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : READ_WRITE_θobs
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      function READ_WRITE_θobs(iSiteName, Path_Date_Output, Path_θ_Input, Path_θ_Output)

         Months = ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]

         # NETCDF
            # Getting daily θ observed
            θdata₀ = Float64.(NetCDF.open(Path_θ_Input, "obsm"))

            N = length(θdata₀)

            θdata = fill(0.0::Float64, N)

            for i=1:N
               θdata[i] =  θdata₀[i] / 100.0
               if θdata[i] < 0.0
                  θdata[i] = NaN
               end
            end

            # Getting the inititilal starting date in the format 1 January 2008   
               Data = NetCDF.open(Path_θ_Input)

               Data2 = NCDatasets.Dataset(Path_θ_Input)

               Date_Start_Obs = Data2["obsm"].attrib["initial_date"]

               # Converting to Year Month Day
                  # Day
                     First_Space = findfirst(" ", Date_Start_Obs)
         
                     Day = Date_Start_Obs[1:First_Space[1]]
                     Date_Start_Obs = replace(Date_Start_Obs, Day => "")
                     Day = parse(Int, Day)

                  # Month
                     Month = 0
                     for iMonth = 1:12
                        if occursin(Months[iMonth],Date_Start_Obs)
                           Month = iMonth
                           Date_Start_Obs = replace(Date_Start_Obs, Months[iMonth] => "")
                           break
                        end
                     end

                  # Year
                     Year = parse(Int, Date_Start_Obs)

                  # Date
                     Date_Start_Obs = Dates.Date(Year,Month,Day)

                # Need to provide dates by adding 1 day
                  Years   = Array{Int64}(undef, N)
                  Months  = Array{Int64}(undef, N)
                  Days    = Array{Int64}(undef, N)
                  Hours   = fill(9::Int64, N)
                  Minutes = fill(0::Int64, N)
                  Seconds = fill(0::Int64, N)

          # CHECKING IF WE HAVE θ IN THE GIVEN DATES        
            # Determening when we start to have θobs with data 
               Add_Date = Date_Start_Obs
               iFirst_NoData = false
               iTrueStart = 1
               for iDay = 1:N
                  if θdata[iDay] < 0 && iFirst_NoData == false
                     iFirst_NoData == true
                     iTrueStart = iDay + 1
                  else
                     iFirst_NoData = true
                  end
                  Years[iDay]  = Dates.year(Add_Date)
                  Months[iDay] = Dates.month(Add_Date)
                  Days[iDay]   = Dates.day(Add_Date)

                  Add_Date     = Date_Start_Obs + Dates.Day(iDay)
               end

               Year_Obs_Start  = Years[iTrueStart]
               Month_Obs_Start = Months[iTrueStart]
               Day_Obs_Start   = Days[iTrueStart]

               Date_Start_Obs = Dates.Date(Year_Obs_Start, Month_Obs_Start, Day_Obs_Start)
                  
            # Determening when we end of having data starting in reverse
               iTrueEnd = N
               iFirst_NoData = false
               for iDay in reverse(1:N)
                  if θdata[iDay] < 0 && iFirst_NoData == false
                     iFirst_NoData == true
                     iTrueEnd = iDay - 1
                  else
                     iFirst_NoData = true
                  end
               end

               Year_Obs_End  = Years[iTrueEnd]
               Month_Obs_End = Months[iTrueEnd]
               Day_Obs_End   = Days[iTrueEnd]

               Date_End_Obs = Dates.Date(Year_Obs_End, Month_Obs_End, Day_Obs_End)

         # SIMULATING DATES      
            # Initial guess
               Year_Sim_Start = 2008 #2008
               Month_Sim_Start = 7 #7
               Day_Sim_Start = 1 #1

               Date_Start_Sim = Dates.Date(Year_Sim_Start, Month_Sim_Start, Day_Sim_Start)

               Year_Sim_End = 2011 #2011
               Month_Sim_End = 12 #12
               Day_Sim_End = 1 #1

               Date_End_Sim = Dates.Date(Year_Sim_End, Month_Sim_End, Day_Sim_End)

            # Checking if the dates are available or else we shift the dates by one year
               if (Date_Start_Sim < Date_Start_Obs)
                  error("$iSiteName Date_Start_Sim $Date_Start_Sim < Date_Start_Obs $Date_Start_Obs")
               end 
                     
               if (Date_End_Sim > Date_End_Obs)
                  Year_Shift = -(Dates.year(Date_End_Sim) - Dates.year(Date_End_Obs) + 1)

                  Year_Sim_Start += Year_Shift
                  Date_Start_Sim = Dates.Date(Year_Sim_Start, Month_Sim_Start, Day_Sim_Start)
                  
                  Year_Sim_End += Year_Shift
                  Date_End_Sim = Dates.Date(Year_Sim_End, Month_Sim_End, Day_Sim_End)

                  if (Date_End_Sim > Date_End_Obs)
                     error("$iSiteName Date_End_Sim $Date_End_Sim > Date_End_Obs $Date_End_Obs")
                  end

                  # Lets check again
                  if (Date_Start_Sim < Date_Start_Obs)
                     error("$iSiteName Date_Start_Sim $Date_Start_Sim < Date_Start_Obs $Date_Start_Obs")
                  end 
               end 

         # CHECKING IF WE HAVE θ on Date_Start_Sim which is used as θᵢₙᵢ 
            # Searching for θᵢₙᵢ
               θᵢₙᵢ = 0.0
               for iDay = 1:N
                  # param.hyPix.Year_Start = 2010
                  if Years[iDay]==Year_Sim_Start && Months[iDay]==Month_Sim_Start && Days[iDay]==Day_Sim_Start
                     θᵢₙᵢ = θdata[iDay]
                  break
                  end
               end

            # Testing
               if θᵢₙᵢ == 0.0
                  error("JulesM error for $iSiteName, date= $(Year_Sim_Start) \\ $(Month_Sim_Start) \\ $(Day_Sim_Start) cannot find θobs")
               end 

               if θᵢₙᵢ < 0.0
                  error("JulesM error for $iSiteName, date= $(Year_Sim_Start) \\ $(Month_Sim_Start) \\ $(Day_Sim_Start) θobs is not available for that date")
                  end 

         # WRITING TO FILES
            # Writting θ
               Header = ["Year";"Month";"Day";"Hour";"Minute";"Second";"Z=200mm"]

               Output = Tables.table( [Years Months Days Hours Minutes Seconds θdata])
         
               CSV.write(Path_θ_Output, Output, header=Header)
            
            # Writting Start end dates
               Header = ["Year_Obs_Start";"Month_Obs_Start";"Day_Obs_Start";"Year_Obs_End";"Month_Obs_End";"Day_Obs_End";"Year_Sim_Start"; "Month_Sim_Start"; "Day_Sim_Start"; "Year_Sim_End"; "Month_Sim_End"; "Day_Sim_End"]

               Output = Tables.table([Year_Obs_Start Month_Obs_Start Day_Obs_Start Year_Obs_End Month_Obs_End Day_Obs_End Year_Sim_Start Month_Sim_Start Day_Sim_Start Year_Sim_End Month_Sim_End Day_Sim_End])
         
               CSV.write(Path_Date_Output, Output, header=Header)


         # println(iSiteName, "=", θᵢₙᵢ)

      return θᵢₙᵢ
      end  # function: READ_WRITE_θobs



   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   #		FUNCTION : READ_WRITE_θobs
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   """
   # Four nc files are saved at each station directory:
   loobos.p2.01d_CP.nc  loobos.p2.01d.nc  loobos.p2.01d_type.nc

   Option" Texture":       loobos.p2.01d_type.nc is the simulation using soil type data from Landcare and Table 2 of Chen and Dudhia (2001) and Brook and Corey model.
   Option" vanGenuchten":  loobos.p2.01d.nc is the simulation using van Genuchten model and the updated soil parameters from Landcare
   Option "BrookCorey":    loobos.p2.01d_CP.nc is the simulation using Brook and Corey model and updated the soil parameters from Landcare.

   For the 4th nc file saved in each station dir such as: sm_sim_obs_11234.nc. It is the daily soil moisture observations and daily soil moisture simulations by JULES using our soil data at Station 11234. Please ignore this simulation for all stations. e.g., sm_sim_obs_11234.nc is for station 11234.

   """
      function READ_WRITE_θJULES(Path_θjules_Input, Path_θjules_Output, Options_θjules)

         OptionsJules = ["Texture","vanGenuchten","BrookCorey"]

         OptionsFile =["loobos.p2.01d_type.nc", "loobos.p2.01d.nc", "loobos.p2.01d_CP.nc"]

         Months = ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]

         # Initializing the dictionary
         i = 1
         OptionsJules_2_OptionsFile = Dict("a"=>"9999") # Initializing
         for iOptionJules in OptionsJules
            OptionsJules_2_OptionsFile[iOptionJules] = OptionsFile[i]
            i += 1
         end

         # Getting θ observed
            Path_θjules = Path_θjules_Input * OptionsJules_2_OptionsFile[Options_θjules]

         # println(NCDatasets.Dataset(Path_θjules))


         θjules = NetCDF.open(Path_θjules, "smcl")

         Time = NetCDF.open(Path_θjules, "time")

         # smcl[x y soil tstep] (1, 1, 4, 4292)

         # println(θjules[1,1,1,1:10])

         Data2 = NCDatasets.Dataset(Path_θjules)

         Date_Start_Obs = Data2[ "timestp"].attrib["time_origin"]
         TimeStep =  Data2[ "timestp"].attrib["tstep_sec"]

      return nothing
      end  # function: READ_WRITE_θobs



      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      #		FUNCTION : READ_JULESSITES
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         function READ_JULES_SITES()
            # Read data
               Data = DelimitedFiles.readdlm(path.JulesMetadata, ',')
            # Read header
               Header = Data[1,1:end]
            # Remove first READ_ROW_SELECT
               Data = Data[2:end,begin:end]
            # Reading
               SiteName, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "SiteName")
         return SiteName
         end  # function: READ_JULESSITES

   end  # module: jules
# ............................................................