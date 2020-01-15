#By Arun & Angel

from ROOT import TFile, TFractionFitter, TObjArray

import pprint

## open root file, created from the SaveHistogramsToRoot step, containing charged hadron isolation distributions which have been grouped into isolated and nonprompt categories
_file = TFile("Isolation_Output.root")

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
    
    ## Add histograms from Isolated and NonPrompt categories to the array 'mc'
    mc = TObjArray(2)
    
    mc.Add(_file.Get("Isolated_"+syst))
    mc.Add(_file.Get("NonPrompt_"+syst))
   
     ## Fit the MC histograms to data 
    fit = TFractionFitter(data,mc,"q")
    
    ## fit.Fit() actually performs the fit
    ## check the fit status
    status = int(fit.Fit())
    
    ## status==0 corresponds to fits that converged, and we can then obtain the fit result 
    fitResults = (fit.GetFitter().Result().Parameters()[0],fit.GetFitter().Result().Parameters()[1])
    
    ## Calculating the scale factor for isolated photons 
    isolatedSF  = data.Integral()*fitResults[0]/mc[0].Integral()
    
    ## Similarly, calculate the scale factor for the nonprompt photons
    nonPromptSF = data.Integral()*fitResults[1]/mc[1].Integral()
    
    ## Calculate the number of events with isolated photons, using the isolatedSF
    isolatedRate = mc[0].GetBinContent(1)*isolatedSF
    ## Calculate the number of events with nonPrompt photons, using the nonPromptSF
    nonPromptRate = mc[1].GetBinContent(1)*nonPromptSF
   
    totalRate = (isolatedRate + nonPromptRate)

    if not status==0:
        print (f"Error in fit while processing {syst} sample: exit status {status}")
    
    ## Now that we know the number of events with isolated photons and the total number of events, we can calculate the photon Purity
    phoPurity = isolatedRate/totalRate

    ## Get the error on the fit parameter for isolated and nonprompt category
    fitError_iso = fit.GetFitter().Result().ParError(0)
    fitError_np = fit.GetFitter().Result().ParError(1)
    
    ## Calculate the error on isolatedRate and nonPromptRate
    isoError = data.Integral()*fitError_iso*mc[0].GetBinContent(1)/mc[0].Integral()
    npError = data.Integral()*fitError_np*mc[1].GetBinContent(1)/mc[1].Integral()

    ## Now we can also calculate the error on photon Purity
    phoPurityErr = ((isoError * (1 + phoPurity) / totalRate)**2 + (npError*phoPurity/totalRate)**2)**0.5

    ## Fill the dictionary "results" with the photonPurity and error in photonPurity for each systematic
    results[syst] = (phoPurity, phoPurityErr)    

    del fit

pp = pprint.PrettyPrinter(indent=4)
pprint.pprint(results)