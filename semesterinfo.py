#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 10:45:35 2021

@author: arest
"""
import sys

default_semester = '2021A'

class semesterinfoclass:
    def __init__(self,semester=None):
        self.setsemester(semester)
        
    def setsemester(self,semester):
        self.semester=semester    

        if semester is None:
            self.programlist=None
            self.program2fieldpattern=None
            return(0)
        
        if semester=='2021A':
            # main programs and assigned night hours
            # YSE to Eta Car: 0.1468+0.6389+1.4861+0.1333+0.6472+1.6667 = 4.71 hours
            # Shen to YSE: 0.6389+0.6472-0.5 = 0.7861 hours
            # DEBASS to YSE: 0.1468+0.1333+0.0535+0.1039+0.0785+0.0344+0.2104+0.1771=0.9379 hours
            YSE2EtaCar = 4.71
            Shen2YSE = 0.7861
            DEBASS2YSE = 0.9379
            Shen2DEBASS = 9.0
            EtaCar2YSE_nextsem = 13.0

            self.programlist = {
                'Shen':7.2+8.233+1.33 -  Shen2YSE - Shen2DEBASS,
                'Martini':1.867,
                'eFEDS':8.75,
                'YSE':70.0 - YSE2EtaCar + Shen2YSE + DEBASS2YSE + EtaCar2YSE_nextsem,
                'DEBASS':29.0 - DEBASS2YSE + Shen2DEBASS,
                'DDF':58.0,
                'DESI':97.0,
                'EtaCar':17.0 + YSE2EtaCar - EtaCar2YSE_nextsem
                }
            
            # patterns to assign Objects 
            self.program2fieldpattern = {
                'Shen':     ['^SN\-C3','^S\-CVZ','SN\-X\d','^CO\d$'],
                'Martini':   ['^E1$','^E3$','^E2$'],
                'eFEDS':     ['^eFEDS'],
                'YSE':       ['^\d\d\d\.\w+\.[abcde]','^YSE$','^\d\d\d\.\w+\.2021'],
                'DEBASS':    ['^2021\w+'],
                'DDF':       ['^COSMOS','^DECaPS.*','^ELAIS'],
                'DESI':      ['^TILEID\:\s+\d+'], 
                'EtaCar':    ['^ec\d\d\d\d'],
                'Drlica_TRADE': ['^DELVE'],
                'Miller_TRADE': ['^n2997'],
                'Rector_TRADE': ['^Cha'],
                'Zenteno_TRADE':['^BLA'],
                'STANDARDS':            ['^E','^SDSS','^LTT','C26202'],
                'TECHSETUP':            ['^pointing','^MaxVis']
                }
            
            self.nights={
                '20210318':1,
                '20210321':1,
                '20210323':1,
#                '20210324':1,
                '20210327':1,
                '20210330':1,
                '20210402':1,
                '20210405':1,
                '20210408':1,
                '20210411':1,
                '20210414':1,
                '20210417':1,
                '20210420':1,
                '20210424':1,
                '20210426':1,
                '20210428':1,
                '20210502':1,
                '20210505':1,
                '20210508':1,
                '20210511':1,
                '20210514':1,
                '20210517':1,
                '20210521':1,
                '20210523':1,
                '20210526':1,
                '20210529':1,
                '20210601':1,
                '20210604':0.5,
                '20210607':0.5,
                '20210610':0.5,
                }
        elif semester=='2021B':
            # main programs and assigned night hours

            self.programlist = {
                'Shen':20.0,
                'Martini':5.0,
                'Liu':20,
                'eFEDS':5.0,
                'YSE':75.0,
                'DEBASS':30.0,
                'DDF':60.0,
                'EtaCar':25.0,
                'Sheppard':5.0
#                'Shen':30.0,
#                'Martini':5.0,
#                'Liu':20,
#                'eFEDS':5.0,
#                'YSE':90.0,
#                'DEBASS':35.0,
#                'DDF':75.0,
#                'EtaCar':25.0,
#                'Sheppard':10.0
                 }
            
            # patterns to assign Objects 
            self.program2fieldpattern = {
                'Shen':      ['^SN\-C3','^S\-CVZ','SN\-X\d','^CO\d$'],
                'Martini':   ['^E1$','^E3$','^E2$'],
                'Liu':       ['BLA'],
                'eFEDS':     ['^eFEDS'],
                'YSE':       ['^\d\d\d\.\w+\.[abcde]','^YSE$','^\d\d\d\.\w+\.2021'],
                'DEBASS':    ['^2021\w+|^2022\w+'],
                'DDF':       ['^COSMOS','^DECaPS.*','^ELAIS'],
                'EtaCar':    ['^ec\d\d\d\d'],
                'Sheppard':  ['twilight'],
                'STANDARDS': ['^E','^SDSS','^LTT','C26202'],
                'TECHSETUP': ['^pointing','^MaxVis']
                }
            
            self.nights={
                '20210916':1.0,
                '20210919':1.0,
                '20210922':1.0,
                '20210925':0.5,
                '20210928':0.5,
                '20211001':0.5,
                '20211004':0.5,
                '20211007':0.5,
                '20211010':0.5,
                '20211013':0.5,
                '20211016':1.0,
                '20211019':1.0,
                '20211022':1.0,
                '20211024':0.5,
                '20211028':0.5,
                '20211031':0.5,
                '20211103':0.5,
                '20211106':0.5,
                '20211109':0.5,
                '20211113':0.5,
                '20211115':0.5,
                '20211118':0.5,
                '20211121':0.5,
                '20211124':0.5,
                '20211127':0.5,
                '20211130':0.5,
                '20211203':0.5,
                '20211206':0.5,
                '20211209':0.5,
                '20211212':0.5,
                '20211215':1.0,
                '20211218':1.0,
                '20211221':1.0,
                '20211223':1.0,
                '20211227':0.5,
                '20211230':0.5,
                '20220102':1.0,
                '20220105':0.5,
                '20220108':0.5,
                '20220112':1.0,
                '20220114':1.0,
                '20220117':1.0,
                '20220120':1.0,
                '20220123':0.5
                }
            sum=0.0
            for k in self.nights:
                sum+=self.nights[k]
            print('total cumulative nights:',sum)
            sys.exit(0)
        elif semester=='2022A':
            # main programs and assigned night hours

            self.programlist = {
                'Shen':20.0,
                'Martini':5.0,
                'Liu':20,
                'eFEDS':5.0,
                'YSE':75.0,
                'DEBASS':30.0,
                'DDF':60.0,
                'EtaCar':25.0,
                'Sheppard':5.0
#                'Shen':30.0,
#                'Martini':5.0,
#                'Liu':20,
#                'eFEDS':5.0,
#                'YSE':90.0,
#                'DEBASS':35.0,
#                'DDF':75.0,
#                'EtaCar':25.0,
#                'Sheppard':10.0
                 }
            
            # patterns to assign Objects 
            self.program2fieldpattern = {
                'Shen':      ['^SN\-C3','^S\-CVZ','SN\-X\d','^CO\d$'],
                'Martini':   ['^E1$','^E3$','^E2$'],
                'Liu':       ['BLA'],
                'eFEDS':     ['^eFEDS'],
                'YSE':       ['^\d\d\d\.\w+\.[abcde]','^YSE$','^\d\d\d\.\w+\.2021'],
                'DEBASS':    ['^2021\w+|^2022\w+'],
                'DDF':       ['^COSMOS','^DECaPS.*','^ELAIS'],
                'EtaCar':    ['^ec\d\d\d\d'],
                'Sheppard':  ['twilight'],
                'STANDARDS': ['^E','^SDSS','^LTT','C26202'],
                'TECHSETUP': ['^pointing','^MaxVis']
                }
            
            self.nights={
                '20220417':1.0,
                '20220420':1.0,
                '20220423':1.0,
                '20220426':1.0,
                '20220429':1.0,
                '20220502':1.0,
                '20220505':1.0,
                '20220508':1.0,
                '20220511':1.0,
                '20220514':1.0,
                '20220517':1.0,
                '20220520':1.0,
                '20220523':1.0,
                '20220526':0.5,
                '20220529':1.0,
                '20220601':0.5,
                '20220604':0.5,
                '20220607':1.0,
                '20220621':1.0,
                '20220624':1.0,
                '20220627':1.0,
                '20220630':1.0,
                '20220703':1.0,
                '20220709':1.0,
                '20220712':1.0,
                '20220715':1.0,
                '20220718':1.0,
                '20220721':1.0,
                '20220724':1.0,
                '20220727':1.0,
                '20220730':1.0,
                 }
            sum=0.0
            for k in self.nights:
                sum+=self.nights[k]
            print('total cumulative nights:',sum)
            #sys.exit(0)
        else:
            raise RuntimeError('%s IS NOT A VALID SEMESTER! ' % semester)
        return(0)
    
