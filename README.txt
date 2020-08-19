READ ME for github repository containing data and code for Demes et al 2020
This repository (which is best opened as an R project) contains two files:

1) Desiccation_data_Aug2020_R2.csv - this is the file containing all of the data for the paper; fields:
	ID - a unique ID assigned to each sample
	Species - macrophyte species from which the sample came
	Clade - clade of species (see Fig 1 in article)
	Habitat - "subtidal" or "intertidal" for macrophytes; "terrestrial" for terrestrial plants
	Drying.time - for desiccation profile tests only: amount of time out of water experienced by the sample 
	RWC - for desiccation profile tests only: relative water content
	Desiccated.Weight.g - for desiccation profile tests only: desiccated weight (following "drying time"; in grams) 
	Tissue.Wet.Weight.g - for desiccation profile tests only: initial wet weight of sample (in grams)
	Water.Lost.g - Initial wet weight of sample minus wet weight post treatment (for desiccation profile only)
	Dry.Weight.g - dry weight of sample
	Treatment - treatment from article (wet [100% RWC], partial [i.e., 75% RWC], dry [0% RWC] OR "Desiccation.profile" for samples used in desiccation profile tests
	Breaking.Force - raw breaking force of working section
	Extensibility - extensibility of sample (measured as proportion of initial length)
	Strength - breaking stress of sample
	Modulus - modulus of elasticty 
	Cross.sectional.area - cross-sectional area of working section where broken (in m^2)
	Thickness_mm - thickness of working section in mm
	Width - width of working section of sample in mm
	Force_per_width - breaking force from article (in Nmm^-1)
	
2) Desiccation_Aug2020_R2.R - this is the R script containing the code used to conduct analyses in the paper
