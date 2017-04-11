#! /usr/bin/env python
# __author__ = 'bachmair'
import argparse
import ConfigParser
import diamond_radiation
import json
import ROOT
import time
import numpy as np
import drawing_class
import root_style
from numpy import arange,linspace
print 'FILE:',__file__
parser = argparse.ArgumentParser(description='Diamond Radiation Hardness calculation.')
parser.add_argument('configfile', metavar='C', nargs='+',
                   help='config file for radiation hardness')
parser.add_argument('-b','--batchmode', action='store_true',
                   help='enables batch mode')
parser.add_argument('-c','--config',help='Main Config File',default = 'config/analysis.cfg')
args = parser.parse_args()


print args.configfile
style = root_style.root_style()
# style.set_style(1200,1200,1)
style.set_style(825,825,1)
main_config = ConfigParser.ConfigParser()
print args.config
main_config.read(args.config)
style.SetBatchMode(args.batchmode)
vary_mfp_ratio = False
if main_config.has_option("Main","vary_mfp_ratio"):
    vary_mfp_ratio = main_config.getboolean("Main","vary_mfp_ratio")
    mfp_range = (0.5,2.44)
    mfp_steps = 10
    if main_config.has_option("Main","mfp_range"):
        mfp_range = json.loads( main_config.get("Main","mfp_range"))
    if main_config.has_option("Main","mfp_steps"):
        mfp_steps = main_config.getint("Main","mfp_steps")
    print mfp_range,mfp_steps
    print  linspace(mfp_range[0],mfp_range[1],mfp_steps)

results = '\n\nRESULTS: \n'
if main_config.has_option("Main","mfp_ratio"):
    print 'MFP ratio: ',main_config.get("Main","mfp_ratio")
for configfile in args.configfile:
    style.postfix = '_ratio_'+main_config.get("Main","mfp_ratio")
    s = '*************** Analyse %s DATA *******************'%configfile
    print '\n'+'*'*len(s)+'\n'+s+'\n'+'*'*len(s)
    diamond_config = ConfigParser.ConfigParser()
    diamond_config.read(configfile)
    diamonds = []
    dias = diamond_config.get('Main', 'diamonds')
    style.prefix = diamond_config.get('Main', 'Directory') + '/'
    try:
        dias = json.loads(dias)
    except Exception as e:
        print "couldn't convert:",dias
        raise e
    print 'Diamonds to be analysed:',dias
    i = 0
    x_max = 0
    x_min = 0
    y_min = 0
    y_max = 0
    for dia in dias:
        items = diamond_config.items(dia)
        # print dia,dict(diamond_config.items(dia))
        d = diamond_radiation.diamond_radiation((dia,dict(diamond_config.items(dia))), main_config, color = style.get_next_color())
        d.analyse()
        diamonds.append(d)

        i+=1

    title = diamond_config.get('Main', 'Title')
    draw_object = drawing_class.drawing_class(style,name=title.replace(' ',''),config=main_config)
    try:
        if True:
            draw_object.Draw_MFP_to_CCD_functions(diamonds)
    except Exception as e:
        print 'error in Draw_MFP_to_CCD_functions',e
        pass
    try:
        pass
        # draw_object.Draw_CCD_to_MFP_conversions(diamonds)
    except Exception as e:
        print 'error in Draw_CCD_to_MFP_conversions',e
        pass
    try:
        pass
        draw_object.Draw_Linerization(diamonds,title)
    except Exception as e:
        print 'error in Draw_Linerization',e
        pass
    # try:
    if True:
        draw_object.DrawAllDamageConstantFit(diamonds,title)
    try:
        draw_object.DrawDamageCurve(diamonds,title)
    except Exception as e:
        print 'error in DrawDamageCurve',e
        raise e


    results += '\t* {title:15s}: {p0:5.3f}({e0:.3f}) x 10^-18 mum^-1 cm^-2, shift: {shift:5.2f} +/- {eshift:5.2f} * 10^15 cm^-2\n'.format(
            title = title,
            p0 = draw_object.get_damage_constant()[0]*1e3,
            e0 = draw_object.get_damage_constant()[1]*1e3,
            shift = draw_object.get_shift(),
            eshift = draw_object.offset_err)
    if vary_mfp_ratio:
        ratio_results = {}
        ratio_results[diamonds[0].mfp_ratio] = [x*1e3 for x in draw_object.get_damage_constant()]
        for mfp_ratio in linspace(mfp_range[0],mfp_range[1],mfp_steps):
            print 'vary mfp range', mfp_ratio
            for d in diamonds:
                d.default_mfp_ratio = mfp_ratio
                d.set_mfp_ratio(mfp_ratio)
                d.analyse()
            ratio_results[mfp_ratio] = draw_object. DrawDamageConstantFit(diamonds,title,plot=False)
            print 'RESULTS of MFP RATIO: %s - %s'%(mfp_ratio,ratio_results[mfp_ratio])
            print '{ratio:4.2f}: {res[0]:.4f} +/- {res[1]:.4f}'.format(ratio=mfp_ratio,res=ratio_results[mfp_ratio])
        draw_object.draw_ratio_variation(ratio_results)
        ratio_items = sorted(ratio_results.items(), key=lambda tup: tup[1][0])
        print ratio_items[0],ratio_items[-1]

        for mfp,res in sorted(ratio_results.items()):
            print '{ratio:4.2f}: {res[0]:.4f} +/- {res[1]:.4f}'.format(ratio=mfp,res=res)

        print ratio_results




print results
