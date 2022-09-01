import Pkg

function PACKAGES(Option_PackageUpdate)
    """This would add the packages automatically if not available"""

    function PACKAGE_MANAGER(Package)
    """ONLY ADD PACKAGE IF NOT AVAILABLE"""
        if  haskey(Pkg.installed(), Package) == false
            println("Adding $Package package because not available...")
            Pkg.add(Package)
            println("$Package package added")
        end
    end # PACKAGE_MANAGER

	# PACKAGES
		PACKAGE_MANAGER("Revise")
		PACKAGE_MANAGER("Tables")
		PACKAGE_MANAGER("BlackBoxOptim")
		PACKAGE_MANAGER("DataFrames")
		PACKAGE_MANAGER("Polynomials")
		# PACKAGE_MANAGER("GRUtils")
		PACKAGE_MANAGER("LaTeXStrings")
		PACKAGE_MANAGER("Optim")
		PACKAGE_MANAGER("Plots")
		PACKAGE_MANAGER("QuadGK")
		PACKAGE_MANAGER("SpecialFunctions")
		PACKAGE_MANAGER("Statistics")
		PACKAGE_MANAGER("PGFPlots")
		PACKAGE_MANAGER("PGFPlotsX")
		PACKAGE_MANAGER("ForwardDiff")
		PACKAGE_MANAGER("CairoMakie")
		PACKAGE_MANAGER("Makie")
		PACKAGE_MANAGER("GLMakie")
		PACKAGE_MANAGER("CSV")
 

			if Option_PackageUpdate
				println("Updating metdata...")
				Pkg.update()
			end
end # PACKAGES

PACKAGES(false)