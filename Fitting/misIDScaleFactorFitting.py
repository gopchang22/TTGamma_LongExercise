from ROOT import TFile, TFractionFitter, TObjArray

import pprint

## open root file, created from the SaveHistogramsToRoot step, containing the e+gamma mass distributions 
_file = TFile("../RootFiles/MisID_Output_electron.root")

## List of systematics
systematics  = ["nominal",
#                "FSRDown",
#                "FSRUp",
#                "ISRDown",
#                "ISRUp",
#                "JERDown",
#                "JERUp",
#                "JESDown",
#                "JESUp",
#                "PDFDown",
#                "PDFUp",
#                "Q2ScaleDown",
#                "Q2ScaleUp",
#                "btagWeight_heavyDown",
#                "btagWeight_heavyUp",
#                "btagWeight_lightDown",
#                "btagWeight_lightUp",
                "eleEffWeightDown",
                "eleEffWeightUp",
                "muEffWeightDown",
                "muEffWeightUp",
#                "puWeightDown",
#                "puWeightUp",
]

results = {}

## Get data from the input root file
data = _file.Get("dataObs")
    
## Loop over the list of systematics    
for syst in systematics:
    
    ## Get histogram from the MisIDele category
    misID = _file.Get("MisIDele_"+syst)
    
    otherMC = _file.Get(f"Other_{syst}")
    ## Add histograms from WGamma and ZGamma categories
    otherMC.Add(_file.Get("WGamma_"+syst))
    otherMC.Add(_file.Get("ZGamma_"+syst))


    mc = TObjArray(2)
    ## Add the histograms from the MC to the mc array
    mc.Add(misID)
    mc.Add(otherMC)

    ## Fit the MC histograms to data 
    fit = TFractionFitter(data, mc, 'q')
    
    ## fit.Fit() actually performs the fit
    ## check the fit status
    status = int(fit.Fit())
    
    ## status==0 corresponds to fits that converged
    
    ## Get the value of fit parameters
    fitResults = fit.GetFitter().Result().Parameters()[0]
    
    ## In order to calculate the electron mis-identification scale factor (SF), we extract the value of the fit parameter for the misID MC and use it to calculate the fraction of mis-identified electrons
    misIDSF  = data.Integral()*fitResults/mc[0].Integral()
    if not status==0:
        print (f"Error in fit while processing {syst} sample: exit status {status}")
        
    ## Fill the dictionary "results" with the misID SF for each systematic
    results[syst] = misIDSF

    del fit

pp = pprint.PrettyPrinter(indent=4)
pprint.pprint(results)
