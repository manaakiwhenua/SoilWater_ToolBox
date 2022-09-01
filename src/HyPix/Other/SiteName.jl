module sitename
   import ..tool
   import DelimitedFiles

   function SITENEAME()

      Path = "D:\\Main\\MODELS\\SoilWater-ToolBox2\\src\\INPUT\\DataHyPix\\JULES\\JULES_LinkingData.csv"

   # = Special for JULES ==
      # Read data
            Data = DelimitedFiles.readdlm(Path, ',')
         # Read header
            Header = Data[1,1:end]
         # Remove first READ_ROW_SELECT
            Data = Data[2:end,begin:end]
         # Reading
            SiteName, ~   = tool.readWrite.READ_HEADER_FAST(Data, Header, "SiteName")

        return SiteName    
   end
end