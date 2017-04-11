import root_style
import ROOT
style = root_style.root_style()
style.set_style(825,825,1)
ROOT.gROOT.SetBatch(0)
xmax = 50e3
thickness = 500
ratio = 1.47
func_string = "y * (1 - y * (1- TMath::Exp(-1/(y))))"

func_string1 = func_string.replace("y", "mfp/(1+[1])")
func_string1 = func_string1.replace("mfp", "(x/[0])")
ccd_e = ROOT.TF1("f_E_to_CCD", func_string1, 0, xmax)

# Holes
func_string2 = func_string.replace("y", "[1]*mfp/(1+[1])")
func_string2 = func_string2.replace("mfp", "(x/[0])")
ccd_h = ROOT.TF1("f_H_to_CCD", func_string1, 0, xmax)

# Electrons + Holes
func_string = func_string1 + "+" + func_string2
func_string = "1/[0]*(%s)" % func_string

ccd_eh = ROOT.TF1("f_mfp_to_ccd_2", func_string, 0, xmax)
ccd_eh.SetParName(0, "Thickness")
ccd_eh.SetParName(1, "#lambda_h/#lambda_e")

ccd_eh.SetParName(0, "Thickness")
ccd_eh.SetParName(1, "#lambda_h/#lambda_e")

ccd_eh_low = ROOT.TF1("f_mfp_to_ccd_low", func_string, 0, xmax)
ccd_eh_up = ROOT.TF1("f_mfp_to_ccd_up", func_string, 0, xmax)

for f in [ccd_e, ccd_h, ccd_eh]:
    if thickness > 0:
        f.SetParameter(0, 1)
    f.SetParameter(1, ratio)

ccd_eh_up.SetLineStyle(2)
ccd_eh_low.SetLineStyle(2)
ccd_eh.SetRange(0,10000)
ccd_eh.SetNpx(100000)

canvas = style.get_canvas('c')
print ccd_e.GetTitle()
print ccd_h.GetTitle()
print ccd_eh.GetTitle()
ccd_eh.Draw()
style.save_canvas(canvas,"mfp_tp_ccd_conversion")