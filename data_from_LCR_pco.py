import datetime
import requests
import json
import numpy as np
import math
import csv
import re

class MET(object):
    """ convert time in MET to a datetime object"""
    mission_start = datetime.datetime(2001,1,1)
    mjd_ref = 51910+7.428703703703703e-4
    def __init__(self, met):
        if met>362793601: met=met-1 # 2012 leap second
        if met>252460801: met=met-1 # 2008 leap second
        if met>157766400: met=met-1 # 2005 leap second
        self.time = MET.mission_start + datetime.timedelta(0,met)
    def __str__(self):
        return str(self.time)


def data_from_LCR(nomesorgente=None,Cadenza=None,TipoFlusso=None,TipoIndice=None,TSminimo=None):
    
    if Cadenza==None:
        Cadenza='monthly' #daily (3d), weekly
    if TipoFlusso==None:
        TipoFlusso='photon' #energy
    if TipoIndice==None:
        TipoIndice='fixed' # ?
    if TSminimo==None:
        TSminimo='4' # 1 2 3 4
    
    if nomesorgente==None:
        nomesorgente4FGL='4FGL%20J1555.7%2B1111'
        nomeassociato='PG 1553+113'
    else:
        nomesorgente4FGL=nomesorgente
        nomesorgente4FGL=nomesorgente4FGL.replace('L J','L%20J')
        nomesorgente4FGL=nomesorgente4FGL.replace('+','%2B')
        nomeassociato=nomesorgente
    
    
    tempoMJD=0
    Flusso=0
    ErroreFlusso=0
    tempoUL=0
    UpperLimitFlusso=0
    maxFlusso=0
    PercentualeUL=0
    
    url ='https'+'://fermi.gsfc.nasa.gov/ssc/data/access/lat/LightCurveRepository/queryDB.php?typeOfRequest=lightCurveData&source_name='+nomesorgente4FGL+'&cadence='+Cadenza+'&flux_type='+TipoFlusso+'&index_type='+TipoIndice+'&ts_min='+TSminimo
    
    Rweb=requests.get(url,allow_redirects=True)
    dati_grezzi = json.loads(Rweb.content)
    
    ################        
    # associate data:  ts flux   flux_upper_limits   flux_error         
    Met_TS=dati_grezzi['ts']
    Met_flux=dati_grezzi['flux']
    Met_error=dati_grezzi['flux_error']
    Met_UL=dati_grezzi['flux_upper_limits']
    
    tempoMJD=np.ones(len(Met_flux))
    Flusso=np.ones(len(Met_flux))
    ErroreFlusso=np.ones(len(Met_flux))
    ErroriRelativi=np.ones(len(Met_flux))
    
    tempoUL=np.ones(len(Met_UL))
    UpperLimitFlusso=np.ones(len(Met_UL))
    
    tempoTS=np.ones(len(Met_TS))
    listaTS=np.ones(len(Met_TS))
    
    mission_start=datetime.datetime(2001,1,1)
    mjd_ref=51910+7.428703703703703e-4
    for i in range(0,len(Met_flux)):
        erroretemporaneo=Met_flux[i][1]-Met_error[i][1]
        if erroretemporaneo>0:
            tempoMJD[i]=(float(Met_flux[i][0])/86400+MET.mjd_ref)
            Flusso[i]=Met_flux[i][1]
            ErroreFlusso[i]=Met_flux[i][1]-Met_error[i][1]
            ErroriRelativi[i]=erroretemporaneo/Met_flux[i][1]
        else:
            tempoMJD[i]=(float(Met_flux[i][0])/86400+MET.mjd_ref)
            Flusso[i]=Met_flux[i][1]
            ErroreFlusso[i]=0
            erroretemporaneo=Met_flux[i][1]-Met_error[i][1]
            ErroriRelativi[i]=erroretemporaneo/Met_flux[i][1]
    
    for i in range(0,len(Met_UL)):
        tempoUL[i]=(float(Met_UL[i][0])/86400+MET.mjd_ref)
        UpperLimitFlusso[i]=Met_UL[i][1]
        
    for i in range(0,len(Met_TS)):
        tempoTS[i]=(float(Met_TS[i][0])/86400+MET.mjd_ref)
        listaTS[i]=Met_TS[i][1]
    
    
    return Flusso,ErroreFlusso,tempoMJD,UpperLimitFlusso,tempoUL