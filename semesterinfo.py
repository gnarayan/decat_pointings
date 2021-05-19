#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 10:45:35 2021

@author: arest
"""

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
            EtaCar2YSE_nextsem = 12.0

            self.programlist = {
                'Shen':7.2+8.233+1.33 -  Shen2YSE,
                'Martini':1.867,
                'eFEDS':8.75,
                'YSE':70.0 - YSE2EtaCar + Shen2YSE + DEBASS2YSE + EtaCar2YSE_nextsem,
                'DEBASS':29.0 - DEBASS2YSE,
                'DDF':58.0,
                'DESI':97.0,
                'EtaCar':17.0 + YSE2EtaCar - EtaCar2YSE_nextsem

#                '2019A-0065_Shen':7.2+8.233+1.33,
#                '2019B-0304_Martini':1.867,
#                '2020A-0906_eFEDS':8.75,
#                '2021A-0037_Shen2':,
#                '2021A-0275_YSE':70.0,
#                '2020B-0053_DEBASS':29.0,
#                '2021A-0113_DDF':58.0,
#                '2021A-0148_DESI':97.0,
#                '2020A-0415_EtaCar':17.0
                }
            
            # patterns to assign Objects 
            self.program2fieldpattern = {
                'Shen':     ['^SN\-C3','^S\-CVZ','SN\-X\d','^CO\d$'],
                'Martini':   ['E1','E3','E2'],
                'eFEDS':     ['^eFEDS'],
                'YSE':       ['^\d\d\d\.\w+\.[abcde]','^YSE$','^\d\d\d\.\w+\.2021'],
                'DEBASS':    ['^2021\w+'],
                'DDF':       ['^COSMOS','^DECaPS.*'],
                'DESI':      ['^TILEID\:\s+\d+'], 
                'EtaCar':    ['^ec\d\d\d\d'],
                'Drlica_TRADE': ['^DELVE'],
                'Miller_TRADE': ['^n2997'],
                'Rector_TRADE': ['^Cha'],
                'Zenteno_TRADE':['^BLA'],
                'STANDARDS':            ['^E','^SDSS','^LTT','C26202'],
                'TECHSETUP':            ['^pointing','^MaxVis']
#                '2019A-0065_Shen':     ['^SN\-C3','^S\-CVZ','SN\-X\d','^CO\d$'],
#                '2019B-0304_Martini':   ['E1','E3'],
#                '2020A-0906_eFEDS':     ['^eFEDS'],
#                '2021A-0037_Shen2':     ['^CO\d$'],
#                '2021A-0275_YSE':       ['^\d\d\d\.\w+\.[abcde]'],
#                '2020B-0053_DEBASS':    ['^2021\w+'],
#                '2021A-0113_DDF':       ['^COSMOS','^DECaPS.*'],
#                '2021A-0148_DESI':      ['^TILEID\:\s+\d+'], 
#                '2020A-0415_EtaCar':    ['^ec\d\d\d\d'],
#                '2019A-0305_Drlica_TRADE': ['^DELVE'],
#                '2021A-0244_Miller_TRADE': ['^n2997'],
#                '2021A-0010_Rector_TRADE': ['^Cha'],
#                '2021A-0149_Zenteno_TRADE':['^BLA'],
#                'STANDARDS':            ['^E','^SDSS','^LTT','C26202'],
#                'TECHSETUP':            ['^pointing','^MaxVis']
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
        else:
            raise RuntimeError('%s IS NOT A VALID SEMESTER! ' % semester)
        return(0)
    