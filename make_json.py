import numpy as np

def individual(json_outpath,json_prefixs,pointRAs,pointDECs,obss,propids,objects,exptypes,programs):
    for json_prefix,pointRA,pointDEC,obs,propid,tobject,exptype,program in zip(json_prefixs,pointRAs,pointDECs,obss,propids,objects,exptypes,programs):
        json_out = []
        json_out.append('[')
        for ob in obs.split(','):
            filt = ob[0]
            exptime = ob[1:]
            json_out.append('\t{')
            json_out.append('\t\t"count": 1,')
            json_out.append('\t\t"filter": "%s",'%filt)
            json_out.append('\t\t"exptime": %.1f,'%float(exptime.strip()))
            json_out.append('\t\t"RA": %f,' % round(float(pointRA), 5))
            json_out.append('\t\t"Dec": %f,' % round(float(pointDEC), 5))
            json_out.append('\t\t"object": "%s",'%tobject)
            json_out.append('\t\t"program": "%s",'%program)
            json_out.append('\t\t"expType": "%s",'%exptype)
            json_out.append('\t\t"note": "None",')
            json_out.append('\t\t"comment": "%s",'%json_prefix)
            json_out.append('\t\t"wait": "False",')
            json_out.append('\t\t"propid": "%s",'%propid)
            json_out.append('\t},')
        json_out.append(']')
        
        np.savetxt(json_outpath+'/'+json_prefix+'.json', json_out, fmt='%s')
        
