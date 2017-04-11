__author__ = 'bachmair'

import json
import ROOT
import array
import math
import numpy as np
from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, \
    FileTransferSpeed, FormatLabel, Percentage, \
    ProgressBar, ReverseBar, RotatingMarker, \
    SimpleProgress, Timer, AdaptiveETA

def json_convert(data):
    ret_val = str(data).replace(' ','')
    return json.loads( ret_val)

class diamond_radiation:
    def __init__(self,data,config,color=1,marker_style=21):
        self.config = config
        self.define_ccd_func()
        if config.has_option('Main','mfp_ratio'):
            self.default_mfp_ratio = config.getfloat('Main','mfp_ratio')
        else:
            self.default_mfp_ratio =1.
        self.set_mfp_ratio(self.default_mfp_ratio)

        if config.has_option('Main','mfp_ratios'):
            self.mfp_ratios = json_convert(config.get('Main','mfp_ratios'))
        if not type(self.mfp_ratios) == list:
            self.mfp_ratios = []
        if not self.default_mfp_ratio in self.mfp_ratios:
            self.mfp_ratios.append(self.default_mfp_ratio)
        # print 'MFP RATIOS:',config.has_option('Main','mfp_ratios'),self.mfp_ratios

        if config.has_option('Main','max_rel_thickness'):
            max_rel_thickness = config.getfloat('Main','max_rel_thickness')
        else:
            max_rel_thickness = .97

        self.marker_style=marker_style
        self.name = data[0]
        self.title = data[1]['name']
        self.type = data[1]['type']
        self.thickness = int(data[1]['thickness'])
        self.fluence_val = json_convert(data[1]['fluence'])
        self.fluence_err = json_convert(data[1]['fluence_err'])
        self.ccd_val = json_convert(data[1]['ccd'])
        self.ccd_val_orig = self.ccd_val
        self.ccd_val = map(lambda x: min(x,self.thickness*max_rel_thickness),self.ccd_val)
        self.ccd_err = json_convert(data[1]['ccd_err'])
        calib_spread = False
        if config.has_option('Main','calib_spread'):
            calib_spread = self.config.getboolean('Main','calib_spread')
        slope_err = False
        if config.has_option('Main','slope_err'):
            slope_err = self.config.getboolean('Main','slope_err')

        self.ccd_calib_spread = [0.] * len(self.ccd_err)
        if 'calib_spread' in data[1]:
            self.ccd_calib_spread = str(data[1]['calib_spread']).replace(' ', '')
            self.ccd_calib_spread = json.loads(self.ccd_calib_spread)

        self.ccd_slope_err = [0.] * len(self.ccd_err)
        if 'slope_err' in data[1]:
            self.ccd_slope_err = str(data[1]['slope_err']).replace(' ', '')
            self.ccd_slope_err = json.loads(self.ccd_slope_err)
            # raw_input( self.ccd_slope_err)

        if slope_err or calib_spread:
            self.ccd_err_orig = self.ccd_err
            self.ccd_err = map(lambda x,y,z:math.sqrt(x**2+y**2+z**2),
                               self.ccd_err_orig,
                               self.ccd_calib_spread,
                               self.ccd_slope_err)
        self.ignore = False
        self.color = color

    def analyse(self):
        self.calculate_damage_constant()
        self.Print()

    def define_ccd_func(self):
        xmax = 50e3
        func_string = "y * (1 - y * (1- TMath::Exp(-1/(y))))"

        # Electrons
        func_string1= func_string.replace("y","mfp/(1+[1])")
        func_string1 = func_string1.replace("mfp","(x/[0])")
        self.ccd_e= ROOT.TF1("f_E_to_CCD",func_string1,0,xmax)

        # Holes
        func_string2= func_string.replace("y","[1]*mfp/(1+[1])")
        func_string2 = func_string2.replace("mfp","(x/[0])")
        self.ccd_h = ROOT.TF1("f_H_to_CCD",func_string1,0,xmax)

        # Electrons + Holes
        func_string = func_string1 +"+" +func_string2
        func_string = "[0]*(%s)"%func_string
        func_string += ';MFP / #mum;CCD / #mum'
        self.ccd_eh = ROOT.TF1("f_mfp_to_ccd_2",func_string,0,xmax)
        self.ccd_eh.SetParName(0,"Thickness")
        self.ccd_eh.SetParName(1,"#lambda_h/#lambda_e")

        self.ccd_eh_low = ROOT.TF1("f_mfp_to_ccd_low",func_string,0,xmax)
        self.ccd_eh_up = ROOT.TF1("f_mfp_to_ccd_up",func_string,0,xmax)
        self.set_mfp_ratio(1,thickness=1)
        self.ccd_eh_up.SetLineStyle(2)
        self.ccd_eh_low.SetLineStyle(2)
        self.ccd_eh_normalized = self.ccd_eh.Clone('ccd_eh_normalized')
        self.ccd_eh_normalized.SetRange(0,40)
        self.ccd_eh_normalized.SetParameter(0, 1)
        self.ccd_eh_normalized.SetParameter(1, 1.47)
        self.ccd_eh_normalized.SetNpx(10000)
        # self.check_functions()
        #
        # self.set_mfp_ratio(2)
        # self.check_functions()
        # self.ccd_eh.SetLineStyle(2)

    def set_all_ratios(self):
        ratio_str = '#frac{#lambda_h}{#lambda_e} = {%f}'
        self.set_mfp_ratio(self.default_mfp_ratio)
        # self.ccd_eh.SetTitle()
        # ratio_str%(self.default_mfp_ratio))
        low =min(self.mfp_ratios)
        self.ccd_eh_low.SetParameter(1,low)
        # self.ccd_eh_low.SetTitle(ratio_str%(low))
        up = max(self.mfp_ratios)
        self.ccd_eh_up.SetParameter(1,up)
        # self.ccd_eh_low.SetTitle(ratio_str%(up))
        self.ccd_eh_up.SetParameter(0,self.thickness)
        self.ccd_eh_low.SetParameter(0,self.thickness)


    def set_mfp_ratio(self,ratio=1.,thickness = -1):
        self.mfp_ratio = ratio
        # print 'new mfp_ratio: ',self.mfp_ratio
        for f in [self.ccd_e,self.ccd_h,self.ccd_eh]:
            if thickness > 0: f.SetParameter(0,1)
            f.SetParameter(1,ratio)

    def check_functions(self):
        print 'Check Functions for ratio: %.2f'%self.mfp_ratio
        for mfp in [0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,10,20,30,40,50,60,70,100]:
            f0 = self.ccd_eh.Eval(mfp)
            f1 = self.ccd_eh.Eval(mfp)
            print '%8.4f | %9.4f | %9.4f | %9.4f'%(mfp,f0,f1,(f0-f1)*100)


    def Print(self):
        print 'Diamond %s with a thickness of %5.1f mum:'%(self.name,self.thickness)
        print 'fluence | CCD [mum] | CCD [mum] | MFP [mum] |'
        print '--------+-----------+-----------+-----------+'
        for i in range(len(self.fluence_val)):
            print '%7.1f | %9.1f | %9.1f |  %9.1f |  %9.1f |  %9.1f'%(self.fluence_val[i],
                                          self.ccd_val_orig[i],
                                          self.ccd_val[i],
                                          self.mfps[i],
                                                                      self.mfp_merr[i],self.mfp_perr[i])
        print

    def GetName(self):
        return self.name

    def mfp_error_est_with_gaus(self,histo,center,cLow=ROOT.kRed,cUp = ROOT.kGreen, nsigma =2):
        xmin = histo.GetXaxis().GetXmin()
        xmax = histo.GetXaxis().GetXmax()
        gausLow = ROOT.TF1('gausMpLow','gaus',xmin,center+1)
        gausUp = ROOT.TF1('gausMpUp','gaus',center-1,xmax)
        gausLow.SetLineColor(cLow)
        gausLow.SetParLimits(1,center-.001,center+.001)
        gausLow.SetParameter(1,center)
        gausUp = ROOT.TF1('gausMeanUp','gaus',center-1,xmax)
        gausUp.SetLineColor(cUp)
        gausUp.SetParLimits(1,center-.001,center+.001)
        gausUp.SetParameter(1,center)

        histo.Fit(gausLow,'NQR','',xmin,center+1)
        histo.Fit(gausUp,'QRN','',center-1,xmax)
        sLow = gausLow.GetParameter(2)
        sUp = gausUp.GetParameter(2)
        histo.Fit(gausLow,'Q+' ,'', max(center-sLow*nsigma,0) ,center+1)
        histo.Fit(gausUp,'Q+','', center-1, min(xmax,center+sUp*nsigma) )
        sGausLow = gausLow.GetParameter(2)
        sGausUp = gausUp.GetParameter(2)
        return gausLow, gausUp ,sGausLow,sGausUp

    def mfp_error_est_mean_with_gaus(self,histo,mfp):
        return self.mfp_error_est_with_gaus(histo,mfp)


    def mfp_error_est_mp_with_gaus(self,histo):
        cLow = ROOT.kBlue
        cUp = ROOT.kOrange
        mp = histo.GetBinCenter(histo.GetMaximumBin())
        return self.mfp_error_est_with_gaus(histo,mp,cLow,cUp)

    def mfp_error_est_with_fwhm(self,histo,bin):
        center = histo.GetBinCenter(bin)
        m = histo.GetBinContent(bin)
        firstBin = histo.FindFirstBinAbove(m/2)
        lastBin = histo.FindLastBinAbove(m/2)
        fwhmLow = histo.GetXaxis().GetBinLowEdge(firstBin)
        fwhmUp  = histo.GetXaxis().GetBinUpEdge(lastBin)
        if  not (fwhmLow < center < fwhmUp):
            print 'Center:   ', center,' @ %d'%bin
            print 'entries:  ',m
            print 'FWHM_low: ',fwhmLow,' @ %d'%firstBin
            print 'FWHM_up:  ',fwhmUp,' @ %d'%lastBin
            raise Exception
        sFwhmLow = center-fwhmLow
        sFwhmUp = fwhmUp - center
        return m,fwhmLow,fwhmUp,sFwhmLow,sFwhmUp

    def mfp_error_est_mean_with_fwhm(self,histo,mfp):
        bin = histo.FindBin(mfp)
        return self.mfp_error_est_with_fwhm(histo,bin)

    def mfp_error_est_mp_with_fwhm(self,histo):
        bin = histo.GetMaximumBin()
        return self.mfp_error_est_with_fwhm(histo,bin)


    def mfp_error_est_with_integral(self,histo,mfp):
        integral = histo.GetIntegral()
        ints = []
        for j in range(histo.GetNbinsX() - 1):
            ints.append(integral[j])
        scale = ints[-1]
        ints = map(lambda x:x/scale,ints)
        thr = (1-.6827)/2.
        print thr
        nLow = next(x[0] for x in enumerate(ints) if x[1] > thr)
        sLow = histo.GetBinCenter(nLow+1)
        nUp = next(x[0] for x in enumerate(ints) if x[1] > 1-thr)
        sUp = histo.GetBinCenter(nUp+1)
        nMean = next(x[0] for x in enumerate(ints) if x[1] > 0.5)
        sMean = histo.GetBinCenter(nMean+1)
        print nLow,nUp, sLow,sUp,mfp,sMean
        sIntegralLow = mfp - sLow
        sIntegralUp  = sUp - mfp
        histo.Integral()
        return sIntegralLow,sIntegralUp

    def convert_to_mfp_simple(self,ccd,ccd_err):
        mfp = self.ccd_eh.GetX(ccd)
        mfpl = self.ccd_eh.GetX(ccd-ccd_err) - mfp
        mfph = self.ccd_eh.GetX(ccd+ccd_err) - mfp
        # print ccd,ccd_err,mfp,mfpl, mfph,'==>', mfp+3*mfpl,mfp+3*mfph
        xmin = mfp+4*mfpl
        xmax = mfp+4*mfph
        sMfpLow = -1*mfpl
        sMfpUp  = mfph
        errors = {
            'Simple': (sMfpLow,sMfpUp)
        }
        histo = None,
        str = ""
        c1 = None
        return histo,str,c1,errors


    def convert_to_mfp(self,ccd,ccd_err):
        if self.config.has_option("MFP conversion","complex"):
            complex_mfp_conversion = self.config.getboolean("MFP conversion","complex")
        else:
            complex_mfp_conversion = False
        if not complex_mfp_conversion:
            return self.convert_to_mfp_simple(ccd,ccd_err)
        nEntries = 1e0
        if self.config.has_option("MFP conversion","nEntries"):
            nEntries = self.config.getfloat("MFP conversion","nEntries")
        if not complex_mfp_conversion:
            nEntries = 1e1
        mfp = self.ccd_eh.GetX(ccd)
        mfpl = self.ccd_eh.GetX(ccd-ccd_err) - mfp
        mfph = self.ccd_eh.GetX(ccd+ccd_err) - mfp
        print 'convert_to_mfp',self.name,ccd,ccd_err,mfp,mfpl, mfph,'==>', mfp+3*mfpl,mfp+3*mfph
        xmin = mfp+4*mfpl
        xmax = mfp+4*mfph
        nbins = 50
        if xmax >  30e3:
            xmax = 30e3
            nbins = 200
        xmin = max(0,xmin)

        print "MFP: ",mfp,"[%05.1f,%05.1f]"%(xmin,xmax)
        if not xmin<mfp<xmax:
            print 'ERROR',xmin,mfp,xmax
        name = 'h_%s_ccd_%.0f'%( self.name,ccd)
        title = '%s: CCD to mfp for %.1f um;mfp / #um; number of entries'%( self.name,ccd)
        histo = ROOT.TH1F(name,title,nbins,xmin,xmax)
        name = 'h_%s_InvMfp_ccd_%.0f'%( self.name,ccd)
        title = '%s: CCD to Inv_mfp for %.1f um;1/mfp / #um^{-1}; number of entries'%( self.name,ccd)
        xmin2 = 1e-6
        x = 1e-8
        xbins  = []
        print 'convert_to_mfp BLA'
        while x < 1:
            xbins.append(x)
            x*=1.01
        xbins.append(1)
        xbins2  = np.array(xbins)
        histo2 = ROOT.TH1F(name,title,len(xbins)-1,xbins2)
        pbar = ProgressBar(widgets=[Percentage(), Bar(),ETA()], term_width=80)
        inv_max = 0
        for i in pbar(range(int(nEntries))):
            y = ROOT.gRandom.Gaus(ccd,ccd_err)
            x = self.ccd_eh.GetX(y)
            histo.Fill(x)
            inv_x = 1./x
            inv_max = max(inv_x,inv_max)
            # print i, x, inv_x
            histo2.Fill(inv_x)
        histo2.GetXaxis().SetRangeUser(1e-8,inv_max)
        mp = histo.GetBinCenter(histo.GetMaximumBin())
        gaus1, gaus2,sGausLow,sGausUp = self.mfp_error_est_mean_with_gaus(histo,mfp)
        gaus3, gaus4,sGaus2Low,sGaus2Up = self.mfp_error_est_mp_with_gaus(histo)

        m,fwhmLow,fwhmUp,sFwhmLow,sFwhmUp = self.mfp_error_est_mean_with_fwhm(histo,mfp)
        m,fwhmMpLow,fwhmMpUp,sFwhmMpLow,sFwhmMpUp = self.mfp_error_est_mp_with_fwhm(histo)
        sIntegralLow,sIntegralUp = self.mfp_error_est_with_integral(histo,mfp)
        sMfpLow = -1*mfpl
        sMfpUp  = mfph
        errors = {
            'Simple': (sMfpLow,sMfpUp),
            'FWHM_mean': (sFwhmLow,sFwhmUp),
            'FWHM_mp': (sFwhmMpLow,sFwhmMpUp),
            'Integral': (sIntegralLow,sIntegralUp),
            'Gaus_mean': (sGausLow,sGausUp),
            'Gaus_mp': (sGaus2Low,sGaus2Up)
        }

        str =  ' \\hline \\num{%5.1f \\pm %3.1f} '%(ccd,ccd_err)
        str += '& \\num{%6.1f} '%mfp
        str += '& \\num{%6.1f} '%mp
        str += '& \\asymUnc{}{%5.1f}{%5.1f} '%(sMfpLow,sMfpUp)
        str += '& \\asymUnc{}{%5.1f}{%5.1f} '%(sFwhmLow,sFwhmUp)
        str += '& \\asymUnc{}{%5.1f}{%5.1f} '%(sFwhmMpLow,sFwhmMpUp)
        str += '& \\asymUnc{}{%5.1f}{%5.1f} '%(sIntegralLow,sIntegralUp)
        str += '& \\asymUnc{}{%5.1f}{%5.1f} '%(sGausLow,sGausUp)
        str += '& \\asymUnc{}{%5.1f}{%5.1f} '%(sGaus2Low,sGaus2Up)
        str += '\\\\'
        print str
        c1 = ROOT.TCanvas()
        histo.SetStats(False)
        # histo.Draw()
        cut = ROOT.TCutG(name+'_mfp_mean_pos',2)
        cut.SetPoint(0,mfp,-1e9)
        cut.SetPoint(1,mfp,+1e9)
        cut.SetLineColor(ROOT.kRed)
        cut.SetLineStyle(2)
        # cut.Draw('same')
        errors['cutMean'] = cut

        cutMp = ROOT.TCutG(name+'_mfp_mp_pos',2)
        cutMp.SetPoint(0,mp,-1e9)
        cutMp.SetPoint(1,mp,+1e9)
        cutMp.SetLineColor(ROOT.kBlue)
        cutMp.SetLineStyle(2)
        # cutMp.Draw('same')
        errors['cutMp'] = cutMp

        gaus11 = gaus1.Clone()
        gaus11.SetLineStyle(2)
        gaus22 = gaus2.Clone()
        gaus22.SetLineStyle(2)
        gaus33 = gaus3.Clone()
        gaus33.SetLineStyle(2)
        gaus44 = gaus4.Clone()
        gaus44.SetLineStyle(2)
        # gaus11.Draw('SAME')
        # gaus22.Draw('SAME')
        # gaus44.Draw('SAME')
        # gaus33.Draw('SAME')
        errors['Gaus_mean_Low'] = gaus11
        errors['Gaus_mean_Up'] = gaus22
        errors['Gaus_mp_Low'] = gaus33
        errors['Gaus_mp_Up'] = gaus44
        errors['Inverse'] = histo2
        histo2.Draw()
        c1.SetLogx()
        c1.Update()

        return histo,str,c1,errors

    def set_thickness(self):
        self.ccd_eh.SetParameter(0,self.thickness)
        self.ccd_eh.SetParameter(0,self.thickness)

    def calculate_mfps(self):
        # ccd(mfp) = 2 x mfp (1 - mfp/t (1-exp(-t/mfp)))
        #x (1-1/2 x (1-Exp[-(2/x)]))

        self.ccd_eh.SetNpx(100000)

        self.set_thickness()
        title =self.ccd_eh.GetTitle()
        # self.ccd_eh.SetTitle('%s;mfp[#mum];ccd[#mum]'%(title))
        self.set_mfp_ratio(self.default_mfp_ratio)
        self.mfps = [self.ccd_eh.GetX(ccd) for ccd in self.ccd_val]
        self.set_mfp_ratio(2)
        self.mfps_up = [self.ccd_eh.GetX(ccd) for ccd in self.ccd_val]
        self.set_mfp_ratio(1)
        self.mfps_low = [self.ccd_eh.GetX(ccd) for ccd in self.ccd_val]

        for i in range(len(self.ccd_val)):
            ratio_low = (self.mfps[i]- self.mfps_low[i])/self.mfps[i]*100.
            ratio_up  = (self.mfps[i]- self.mfps_up[i] )/self.mfps[i]*100.
            # print '%8.2f: %8.2f, %8.2f, %8.2f'%(self.ccd_val_orig[i],self.mfps[i], ratio_low,ratio_up)
        if self.default_mfp_ratio == 1:
            print ratio_low, ratio_up
        self.inv_mfps = [1/x for x in self.mfps]
        self.mfp_convs =[]
        for ccd, ccd_err in zip(self.ccd_val, self.ccd_err):
                self.mfp_convs.append( self.convert_to_mfp(ccd,ccd_err))
        mccd = map(lambda x,y: x-y,self.ccd_val,self.ccd_err)
        self.mfp_merr = map(lambda x,y:abs(x-y),[self.ccd_eh.GetX(ccd) for ccd in mccd],self.mfps)
        self.inv_mfp_merr = map(lambda x,y: abs(1/x**2*y),self.mfps,self.mfp_merr)

        pccd = map(lambda x,y: x+y,self.ccd_val,self.ccd_err)
        self.mfp_perr = map(lambda x,y:abs(x-y),[self.ccd_eh.GetX(ccd) for ccd in pccd],self.mfps)
        self.inv_mfp_perr= map(lambda x,y: abs(1/x**2*y),self.mfps,self.mfp_perr)


        #devide by to to get mfp_e/mfp_h
        #self.inv_mfps = map(lambda x: x/2, self.inv_mfps)
        #self.inv_mfp_merr = map(lambda x: x/2, self.inv_mfp_merr)
        #self.inv_mfp_perr = map(lambda x: x/2, self.inv_mfp_perr)
        # print 'mfps of %s: '%self.name
        # for i in range(len(self.mfps)):
        #     print '\t%4.1f+/- %3.1f: %6.2f + %6.2f - %6.2f'%(self.fluence_val[i],self.fluence_err[i],self.mfps[i],self.mfp_perr[i],self.mfp_merr[i]),
        #     print ' %6.2e + %6.2e - %6.2e'%(self.inv_mfps[i],self.inv_mfp_merr[i],self.inv_mfp_perr[i])
        self.ccd_eh.SetRange(0,10000)

    def validate_values(self):
        values = [self.fluence_val,self.fluence_err,self.ccd_val,self.ccd_err]
        le = map(lambda x:len(x),values)
        s = set(le)
        l = len(s)
        if  l != 1:
            print 'Values:',values
            print 'Lengths:',le
            print 'Set: ',s,l
            raise Exception('Length of data vectors do not agree for sample "%s"'%self.name)
        return True

    def get_max_fluence(self):
        return max(map(lambda x,y: x+y, self.fluence_val,self.fluence_err))

    def get_min_fluence(self):
        return min(map(lambda x,y: x-y, self.fluence_val,self.fluence_err))

    def get_min_inv_mfp(self):
        return max(map(lambda x,y: x-y, self.inv_mfps,self.inv_mfp_merr))

    def get_max_inv_mfp(self):
        return max(map(lambda x,y: x+y, self.inv_mfps,self.inv_mfp_perr))


    def calculate_damage_constant(self):
        self.validate_values()
        self.calculate_mfps()
        self.gr = ROOT.TGraphAsymmErrors(len(self.fluence_val),
                                    np.array(self.fluence_val),np.array(self.inv_mfps),
                                    np.array(self.fluence_err),np.array(self.fluence_err),
                                    np.array(self.inv_mfp_merr),np.array(self.inv_mfp_perr))
        self.gr.SetName(self.name)
        self.gr.SetTitle(self.name)
        self.gr.GetXaxis().SetTitle('Fluenece [10^{15} cm^{-2}]')
        self.gr.GetYaxis().SetTitle('1/mfp [1/#mum]')
        self.gr.GetYaxis().SetTitleOffset(1.2)
        delta = max(self.mfps) - min(self.mfps)
        self.fit = ROOT.TF1("fit","pol1",min(self.mfps)-.1*delta,20)#max(self.mfps)+.1*delta)
        self.fit.SetLineColor(self.color)
        self.fit.SetLineStyle(2)
        self.fit.SetFillColor(0)
        self.gr.SetLineColor(ROOT.kBlack)
        self.gr.SetMarkerColor(self.color)
        self.gr.SetMarkerStyle(self.marker_style)
        self.gr.SetFillColor(0)
        self.gr.Fit(self.fit,'Q')
        # print 'Fit of {name}: {p0:f} +/- {e0:f}   | {p1:f} +/- {e1:f} '.format(
        #     name = self.name,
        #     p0 = self.fit.GetParameter(0),
        #     e0 = self.fit.GetParError(0),
        #     p1 = self.fit.GetParameter(1),
        #     e1 = self.fit.GetParError(1)
        # )
        self.print_fit_results()
        # print '{name} & \SI{{ {p0:f} \pm {e0:f} }} {{ X }}}}& \SI{{ {p1:f} \pm {e1:f} }} {{ X }} \\'\

        self.fit_long = self.fit.Clone('fit_clone')
        self.fit_long.SetRange(0,30)
        self.gr.GetListOfFunctions().Add(self.fit_long)
        retVal = ((self.fit.GetParameter(0),self.fit.GetParError(0)),
                  (self.fit.GetParameter(1),self.fit.GetParError(1)),
                  (self.fit.GetChisquare(),self.fit.GetNDF()))
        self.fit_val = retVal
        if self.default_mfp_ratio == 1:
            # print retVal
            pass
        return retVal
        #calculate
        pass

    def print_fit_results(self):
        print '\t{name} & {p0:.2f} \pm {e0:.2f} & {p1:.2f} \pm {e1:.2f} \\\\'.format(
            name = self.title,
            p0 = self.fit.GetParameter(0)*1e3,
            e0 = self.fit.GetParError(0)*1e3,
            p1 = self.fit.GetParameter(1)*1e3,
            e1 = self.fit.GetParError(1)*1e3
        )