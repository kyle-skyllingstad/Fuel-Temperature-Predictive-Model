# Fuel-Temperature-Predictive-Model

# Overview
This predictive framework calculates temperatures under which chemical mixtures combust most efficiently. This program is hardcoded to predict this optimal combustion temperature, or superheat limit, for a particular surrogate (replacement mixture) for a given jet fuel. Its hardcoded inputs are commonly-known characteristics of each of the different components in the mixture, and it uses predictive theory to calculate properties of the mixture including its superheat limit. This framework is included in the publication, "Improving performance of multi-hole fuel injectors by predicting flash boiling conditions for Jet-A and its surrogate" published in Energy & Fuels, 2019. 

# How to Use
This predictive model consists of three separate programs that together calculate the superheat limit of the chosen surrogate. They are: PlotNucTempPseudo.m, get3propertiesnew3.m, and MixtureCriticalTP.m. To run the package, simply run the main host program, PlotNucTempPseudo.m. No data files are necessary to run this bundle. It should generate a plot that illustrates the distribution of surrogate superheat limits across a variety of nucleation rates (essentially fuel heating rates). 
