README
Title: “Does the @realDonaldTrump really matter to financial markets?”

Authors: Allyson L. Benton and Andrew Q. Philips.

Journal: American Journal of Political Science. 

Last updated: 9/17/19

Description: replication code, data, and codebook needed to reproduce all results in the main manuscript and supplemental materials

Notes: 
	-Please note that due to differences in the operating system and types of Stata, results will differ very slightly. This appears to cause different results for those reported in Table 2, Model 2 and Table 3, Model 10, although the log likelihoods remain the same. Users trying to get the exact same output should use Stata SE (version 15.1) for Mac
	—users will need to set their Stata graphics scheme to “burd” to replicate figures (via ”findit scheme burd”)
	—users will need to have user-written “kpss” downloaded to replicate Table 1 in the SI (“findit kpss”)
	—users will need to have user-written “sutex” downloaded to replicate Table 2 in the SI (“findit sutex”)
	—files use Stata’s ‘preserve’, ‘restore’ and ‘local’ functions for creating simulation figures. Users will need to run all lines as a single do within a figure to recreate it (e.g., to run Figure 3, highlight and run lines 115-219
	—the x-axis labels in Figure 2b were done by hand using Stata’s graph editor

Files:
	—replication-main.do: code to replicate all tables and figures in the manuscript
	—replication-SI.do: code to replicate all tables and figures in the supplemental materials
	—tweetpeso-final.dta: Stata dataset needed to run “replication-main.do” and “replication-SI.do”. Can be read by Stata 11+ (though all analyses were originally conducted using version 15.1)
	—codebook.pdf: codebook for “tweetpeso-final.dta”
	