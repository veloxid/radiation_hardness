__author__ = 'bachmair'

import copy
import ROOT
import math
import time
import numpy as np
import root_style
import diamond_radiation

class drawing_class:
    def __init__(self,style,name=''):
        self.style= style
        self.canvas = style.get_canvas('Linerization'+name)
        self.offset_phi = 0
        self.prefix = self.style.prefix
        self.postfix = ""
        self.fitsDamageConstant = {}
        self.grsDamageConstant = {}
        # print 'PREFIX: ',self.prefix
        pass

    def Draw_CCD_to_MFP_conversions(self,diamonds):
        print 'Draw_CCD_to_MFP_conversions'
        self.canvas.SetTitle('CCD to MFP conversion')
        self.canvas.cd()
        prefix = self.style.prefix
        try:
            for d in diamonds:
                d_name = d.name.replace(' ','_')
                self.style.prefix ='%s/%s_'%(d_name,d_name)
                print '',d_name,self.style.prefix,prefix
                for conv in d.mfp_convs:
                    #conv = histo,str,c1,cut,gaus11,gaus22
                    #          0 , 1 , 2, 3 , 4    ,  5
                    print conv
                    errors = conv[-1]
                    h = conv[0]
                    s = conv[1]
                    cut1 = errors['cutMean']  if errors.has_key('cutMean') else None
                    cut2 = errors['cutMp'] if errors.has_key('cutMP') else None
                    g1 = errors['Gaus_mean_Low'] if errors.has_key('Gaus_mean_Low') else None
                    g2 = errors['Gaus_mean_Up'] if errors.has_key('Gaus_mean_Up') else None
                    g3 = errors['Gaus_mp_Low'] if errors.has_key('Gaus_mp_Low') else None
                    g4 = errors['Gaus_mp_Up'] if errors.has_key('Gaus_mp_Up') else None
                    print h
                    h.Draw()
                    if cut1: cut1.Draw('same')
                    if cut2: cut2.Draw('same')
                    if g1: g1.Draw('same')
                    if g2: g2.Draw('same')
                    if g3: g3.Draw('same')
                    if g4: g4.Draw('same')
                    self.canvas.Update()
                    self.style.save_canvas(self.canvas,'CCD_to_MFP_'+h.GetName())

                    print d_name
                    print s
        except Exception as e:
            print 'Error in  Draw_CCD_to_MFP_conversions: ',e
        self.style.prefix=prefix

    def Get_Data_Graph(self,diamond):
        mfp_perr = diamond.mfp_perr
        mfp_merr = diamond.mfp_merr
        mfps = diamond.mfps
        ccd = diamond.ccd_val
        ccd_merr = diamond.ccd_err
        ccd_perr = diamond.ccd_err
        nPoints  = len(ccd)
        g = ROOT.TGraphAsymmErrors(nPoints,np.array(mfps),np.array(ccd),
                                   np.array(mfp_merr),np.array(mfp_perr),
                                   np.array(ccd_merr),np.array(ccd_perr))
        g.SetMarkerStyle(5)
        return g

    def Draw_MFP_to_CCD_functions(self,diamonds):
        self.canvas.SetTitle('MFP to CCD conversion')
        self.canvas.cd()
        try:
            for d in diamonds:
                print 'Draw_MFP_to_CCD_functions',d.name
                d.set_all_ratios()
                name = d.name.replace(' ','_')
                g = self.Get_Data_Graph(d)
                self.style.prefix = self.prefix+name+'/'
                cut = ROOT.TCutG('thickness_%s'%name,2)
                cut.SetPoint(0,-1e10,d.thickness)
                cut.SetPoint(1,1e10,d.thickness)
                cut.SetLineStyle(2)
                cut.SetLineColor(ROOT.kRed)
                d.ccd_eh.Draw()
                d.ccd_eh.GetHistogram().SetMaximum(d.thickness*1.1)
                d.ccd_eh.GetHistogram().GetXaxis().SetTitle('MFP / #mum')
                d.ccd_eh.GetHistogram().GetYaxis().SetTitle('CCD / #mum')
                d.ccd_eh_low.Draw('same')
                d.ccd_eh_up.Draw('same')
                cut.Draw('Same')
                g.Draw('PE')
                self.canvas.Update()
                self.style.save_canvas(self.canvas,'MFP_to_CCD_'+name)

                d.ccd_eh.SetRange(0,5e4)
                d.ccd_eh.Draw()
                d.ccd_eh.GetHistogram().SetMaximum(d.thickness*1.1)
                d.ccd_eh_low.Draw('same')
                d.ccd_eh_up.Draw('same')
                cut.Draw('Same')
                g.Draw('PE')
                self.canvas.SetLogx()
                self.style.save_canvas(self.canvas,'MFP_to_CCD_logx_'+name)
                self.canvas.SetLogx(False)
                d.ccd_eh.SetRange(0,1e3)
                d.ccd_eh.Draw()
                d.ccd_eh.GetHistogram().SetMaximum(d.thickness*1.1)
                d.ccd_eh_low.Draw('same')
                d.ccd_eh_up.Draw('same')
                func = ROOT.TF1("func","pol1",0,1e3)
                func.SetParameter(0,0)
                func.SetParameter(1,1)
                func.SetLineStyle(2)
                func.SetLineColor(ROOT.kBlack)
                func.Draw('same')
                cut.Draw('Same')
                g.Draw('PE')
                self.canvas.Update()
                self.style.save_canvas(self.canvas,'MFP_to_CCD_zoom_'+name)
        except Exception as e:
            print 'Cannot Draw MFP to CCD ',e
        self.style.prefix = self.prefix

    def Draw_Linerization(self,diamonds,title):
        x_max = 0
        x_min = 0
        y_max = 0
        y_min = 0

        ROOT.gStyle.SetOptFit(0)
        for dia in diamonds:
            self.style.prefix = self.prefix+dia.name+"/"
            # dia.gr.SetLineColor(i+1)
            dia.gr.Draw('PA')
            self.canvas.Update()
            self.style.save_canvas(self.canvas,'DamageLinearization_%s'%dia.name)

        self.style.prefix = self.prefix
        for d in diamonds:
            x_min = min(d.get_min_fluence(),x_min)
            x_max = max(d.get_max_fluence(),x_max)
            y_min = min(d.get_min_inv_mfp(),y_min)
            y_max = max(d.get_max_inv_mfp(),y_max)
        i = 0
        frame = self.canvas.DrawFrame(x_min-.05*(x_max-x_min),
                     y_min-.05*(y_max-y_min),
                     x_max+.05*(x_max-x_min),
                     y_max-.05*(y_max-y_min),
                     '%s;Fluence #times 10^{-15} / #mum^{-1}cm^{-2};1/MFP   /  1/#mum'% title)
        frame.GetYaxis().SetTitleOffset(1.9)
        self.canvas.SetTitle('Linerization')

        ROOT.gStyle.SetOptFit(0)
        for dia in diamonds:
            # dia.gr.SetLineColor(i+1)
            dia.gr.Draw('P')

            i+=1
            self.canvas.Update()
            # print '%10s:'%dia.name,'%5.5f +/- %5.5f'%(dia.fit_val[1][0],dia.fit_val[1][1])
        leg = self.style.make_legend(.2,.95,len(diamonds))
        for d in diamonds:
            leg.AddEntry(d.gr,d.title,'LP')
        leg.Draw()
        self.canvas.Update()
        # raw_input()
        self.style.save_canvas(self.canvas,'DamageLinearization')

    def DrawAllDamageConstantFit(self,diamonds,title):
        types = set([dia.type for dia in diamonds])
        for t in types:
            dias = filter(lambda x:x.type == t,diamonds)
            self.DrawDamageConstantFit(dias,title,type = t)
        self.DrawDamageConstantFit(diamonds,title,type = 'all')
        # raw_input()
        return

    def DrawDamageConstantFit(self,dias,title, type = 'all'):
        print 'Draw Damage Constant Fit',type, 'for ',[d.name for d in dias]
        constant = []
        constant_err = []
        self.canvas.Clear()
        self.canvas.Update()
        self.canvas.SetTitle('Damage Constant Fit')
        for dia in dias:
            constant.append(dia.fit_val[1][0]*1e-15/1e-18)
            constant_err.append(dia.fit_val[1][1]*1e-15/1e-18)

        x = range(len(constant))
        x = [float(i+1) for i in x]
        ex = [1e-5 for i in x]
        ROOT.gStyle.SetOptFit(1)
        # print 'constants', len(constant)
        gr = ROOT.TGraphErrors(len(x),np.array(x),np.array(constant),np.array(ex),
                            np.array(constant_err))
        gr.SetName('gDamageConstantFit_{0}'.format(type))
        gr.SetTitle('Damage Constant of %s;;k_{mfp} #times 10^{-18} / #mum^{-1}cm^{-2}'%title)
        gr.Draw('AP')
        gr.GetYaxis().SetTitleOffset(1.9)
        fit = ROOT.TF1('fit_{0}'.format(type),'pol0',min(x),max(x))
        gr.Fit(fit,'Q')
        # print gr.GetN()
        hint = gr.Clone()
        ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(hint)
        if type == 'all':
            self.fitDamageConstant = fit
            self.grDamageConstant  = gr

        self.fitsDamageConstant[type] = fit
        self.grsDamageConstant[type]  = gr
        # hint.SetStats(False)
        hint.SetFillColor(18)
        hint.SetFillStyle(3144)
        print 'Draw',gr.GetName(),'with',gr.GetN()
        gr.Draw('APL')
        gr.Draw('AP')
        hint.Draw("e3 same")
        gr.Draw('PE')
        self.canvas.Update()
        xax = gr.GetXaxis()
        print 'x: ',xax.GetXmin(),xax.GetXmax(),xax.GetNbins()
        for i in range(len(dias)):
            bin_index = xax.FindBin(i+1)
            xax.SetBinLabel(bin_index,dias[i].title)
            print i,bin_index,xax.GetBinCenter(bin_index),dias[i].title, gr.GetX()[i],gr.GetY()[i]
        ROOT.gPad.Modified()
        ROOT.gPad.Update()
        gr.Draw('APL')
        gr.Draw('AP')
        hint.Draw("e3 same")
        gr.Draw('PE')
        self.canvas.Update()
        save_as = 'DamageConstantFit_{type}'.format(type=type)
        print 'save_as',save_as
        self.style.save_canvas(self.canvas,save_as)
        raw_input()
        return fit.GetParameter(0),fit.GetParError(0)

    @staticmethod
    def GetSortedDiamonds(diamonds):
        singles = []
        polys = []
        for dia in diamonds:
            if dia.type == 'single':
                singles.append(dia)
            else:
                polys.append(dia)
        return singles,polys

    @staticmethod
    def get_mean_and_error(values):
        mean = 0
        sigma = 0
        n = 0
        is_tuples = 0
        for val in values:
            if isinstance(val, (list, tuple)):
                if is_tuples == -1:
                    raise Exception('Cannot merge different types')
                is_tuples = 1
            else:
                if is_tuples == 1:
                    raise Exception('Cannot merge different types')
                is_tuples = -1
        for val in values:
            if is_tuples == 1:
                is_tuples = True
                mean += val[0] #/val[1]**2
                sigma += val[0]**2
            else:
                mean += val
                sigma += val**2
            n += 1
        if n ==0:
            mean = 0
            sigma = 0
        elif sigma == 0:
            sigma = 0
        elif is_tuples == 0:
            mean = mean/sigma
            sigma = 1/sigma
            sigma = math.sqrt(sigma)
        else:
            mean/= n
            sigma /= n
            sigma = math.sqrt(sigma - mean**2)
        # print values,'==>',mean,sigma
        if is_tuples == 1 and len(values) == 1:
            sigma = values[0][1]
        return mean, sigma

    def get_damage_constant(self):
        return self.fitDamageConstant.GetParameter(0)/1e3, self.fitDamageConstant.GetParError(0)/1e3

    def ExtractFluenceOffset(self,singles,polys):
        print 'ExtractFluenceOffset'
        c_singles = [x.fit_val[0] for x in singles]
        c_polys = [x.fit_val[0] for x in polys]
        # c_polys_err = [x.fit_err[0] for x in polys]
        verb  = False
        if verb: print 'SINGLES: ',c_singles
        if verb: print 'POLIES:  ',c_polys
        c_single = self.get_mean_and_error(c_singles)
        c_poly = self.get_mean_and_error(c_polys)
        damage_constant, damage_error = self.get_damage_constant()
        if verb: print 'c_single',c_single
        if verb: print 'c_poly  ',c_poly
        if verb: print 'damage  ',damage_constant
        offset_phi = 0
        offset_err = 0
        if damage_constant !=0:
            offset_phi = (c_poly[0]-c_single[0])/damage_constant
            offset_err  = (c_poly[1]/damage_constant)**2
            offset_err += (c_single[1]/damage_constant)**2
            offset_err += ((c_poly[0]-c_single[0])/damage_constant**2 *damage_error)**2
            offset_err = math.sqrt(offset_err)
        if verb: print 'offset   ',offset_phi,offset_err
        if verb: print 'retvals:',offset_phi, c_poly,c_single,damage_constant
        if verb: raw_input()
        self.offset_phi = offset_phi
        self.offset_err = offset_err
        print offset_err,offset_phi
        return offset_phi,c_single

    def get_shift(self):
        return self.offset_phi

    def DrawDamageCurve(self,diamonds,title):
        print 'DrawDamageCurve'

        singles,polys = self.GetSortedDiamonds(diamonds)
        offset_phi,constant = self.ExtractFluenceOffset(singles,polys)

        self.graphs = []
        x_min = 1e10
        x_max = -1e10
        y_min = 1e10
        y_max = -1e10
        self.canvas.SetTitle('Damage Curve')
        self.canvas.cd()
        self.canvas = self.style.get_canvas('Damage Curve')
        for dia in diamonds:
            dia.print_fit_results()
            vx = dia.fluence_val
            if dia.type != 'single':
                vx = map(lambda x: x+offset_phi,vx)
            vy = dia.mfps
            ex = dia.fluence_err
            eyl = dia.mfp_merr
            eyh = dia.mfp_perr
            gr = ROOT.TGraphAsymmErrors(len(vx),
                                    np.array(vx),np.array(vy),
                                    np.array(ex),np.array(ex),
                                     np.array(eyl),np.array(eyh))
            gr.SetName('gDamageCurve_'+dia.name)
            gr.SetTitle(dia.gr.GetTitle())
            gr.SetMarkerColor(dia.color)
            gr.SetMarkerStyle(20)
            self.canvas.cd()
            gr.Draw('APC')
            self.canvas.Update()
            self.graphs.append(copy.deepcopy(gr))
            x_max = max(x_max,max(vx))
            x_min = min(x_min,min(vx))
            y_max = max(y_max,max(vy))#gr.GetHistogram().GetMaximum())
            y_min = min(y_min,min(vy))
        print 'Result &  -- & {p0:.2f} \pm {e0:.2f} \\\\'.format(
            p0 = self.fitDamageConstant.GetParameter(0),
            e0 = self.fitDamageConstant.GetParError(0)
        )
        dx = .05*(x_max-x_min)
        dy = .05*(y_max-y_min)
        y_min = .9 * y_min
        y_min = max(y_min,1)
        self.canvas.Clear()
        self.canvas.Update()
        self.canvas.cd()

        damage_constant, damage_constant_err = self.get_damage_constant()
        print  damage_constant, damage_constant_err
        # damage_constant_err = self.fitDamageConstant.GetParError(0)*2/1e3
        self.fit = ROOT.TF1('damge_curve',"[0]/(1+[1]*[0]*x)",0,20)
        self.fit.SetParName(0,'#lambda_0')
        self.fit.SetParName(1,'k')
        self.fit.FixParameter(1,damage_constant)
        self.fit.SetParError(1,damage_constant_err)
        self.fit.SetLineColor(self.style.get_next_color())
        self.fit.SetNpx(1000)
        self.fit.SetParameter(0,constant[0])
        self.fit_high = self.fit.Clone('fit_high')
        self.fit_high.FixParameter(1,damage_constant+damage_constant_err)
        self.fit_high.SetLineStyle(3)
        self.fit_low = self.fit.Clone('fit_low')
        self.fit_low.FixParameter(1,damage_constant-damage_constant_err)
        self.fit_low.SetLineStyle(3)
        self.fit.SetParameter(0,constant[0])
        # print constant[0],1/constant[0]
        self.fit.SetParError(0,constant[1])
        print "FIT:"
        for i in range(2):
            print '\t* ',i,'%5.2e +/- %5.2e'%(self.fit.GetParameter(i),self.fit.GetParError(i))
        self.cut = ROOT.TCutG('Shift',2)
        self.cut.SetName('Shift')
        self.cut.SetTitle('Shift: %.1f#times 10^{15} cm^{-2}'%(offset_phi))
        self.cut.SetLineColor(ROOT.kRed)
        self.cut.SetLineStyle(2)
        self.cut.SetPoint(0,offset_phi,-1e10)
        self.cut.SetPoint(1,offset_phi,+1e10)

        self.canvas.cd()
        frame_title = title+';Fluence / #times 10^{15} cm^{-2};mfp /#mum'
        frame = self.canvas.DrawFrame(x_min-dx,y_min,x_max+dx,y_max+dy,frame_title)
        self.canvas.Update()
        frame.GetYaxis().SetTitleOffset(1.9)

        print 'fitting:'
        for gr in self.graphs:
            print gr.GetName(),gr.GetTitle()
            frame.Draw()
            gr.Fit(self.fit,'QN')
            print '  fitted: ','%10s'%gr.GetName(),self.fit.GetParameter(0)
            gr.Draw('P')
            self.fit.Draw('same 4')
            self.fit.Draw('SAME')
            self.fit_high.Draw('SAME')
            self.fit_low.Draw('SAME')
            self.canvas.Update()
            try:
                dia = diamonds[self.graphs.index(gr)]
                dia_name = dia.name
                self.style.prefix = self.prefix+dia_name+'/'
                title = 'DamageCurve_'+dia_name
                print 'save',title
                self.style.save_canvas(self.canvas,title)
                self.canvas.SetLogy()
                self.canvas.Update()
                self.style.save_canvas(self.canvas,title+'_logy')
                self.canvas.SetLogy(False)
            except Exception as e:
                print 'Exception:',e
        self.style.prefix = self.prefix
        frame.Draw()
        for gr in self.graphs:
            # gr.Fit(self.fit,'QN')
            print '  fitted: ','%10s'%gr.GetName(),self.fit.GetParameter(0)
            gr.Draw('P')
        self.fit.Draw('same 4')
        self.fit.Draw('SAME')
        self.fit_high.Draw('SAME')
        self.fit_low.Draw('SAME')
        self.cut.Draw('same')
        self.leg = self.style.make_legend(.45,.93,len(diamonds)+2)
        self.leg.AddEntry(self.cut,self.cut.GetTitle(),'L')
        for d in diamonds:
            t = d.title
            if d.type != 'single':
                t += ' shifted'
            self.leg.AddEntry(d.gr,t)
        self.leg.AddEntry(self.fit,'Fit','L')
        self.leg.AddEntry(None,'k_{mfp}: %2.2e #mum^{-1} cm^{-2}'%(damage_constant/1e15),'')
        self.leg.Draw()
        self.canvas.Update()
        self.style.save_canvas(self.canvas,'DamageCurve')
        self.canvas.SetLogy()
        self.canvas.Update()
        self.style.save_canvas(self.canvas,'DamageCurveLog')
        print
