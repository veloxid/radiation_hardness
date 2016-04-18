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
style.set_style(1200,1200,1)
main_config = ConfigParser.ConfigParser()
print args.config
main_config.read(args.config)
style.SetBatchMode(args.batchmode)
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
    draw_object = drawing_class.drawing_class(style,name=title.replace(' ',''))
    try:
        if False:
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
    try:
        draw_object.DrawAllDamageConstantFit(diamonds,title)

        pass
    except Exception as e:
        raise e
        print 'error in DrawDamageConstantFit',e
        pass
    try:
        draw_object.DrawDamageCurve(diamonds,title)
    except Exception as e:
        print 'error in DrawDamageCurve',e
        raise e
    results += '\t* %15s: %5.2f x 10^-18 mum^-1 cm^-2, shift: %5.2f +/- %5.2f * 10^15 cm^-2\n'%(title,
                                                                                                draw_object.get_damage_constant()[0]*1e3,
                                                                                                draw_object.get_shift(),
                                                                                                draw_object.offset_err)
print results
