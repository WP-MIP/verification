#!/usr/bin/env python
import glob, os, sys, json
import earthkit as ek
import xarray as xr
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import rpnpy.librmn.all as rmn
import spharm
    
class Spec(dict):

    def __init__(self, centre, stream, ftype, version, field, level, lead):
        self['centre'] = centre
        self['stream'] = stream
        self['ftype'] = ftype
        self['version'] = "%02d" % int(version)
        self['field'] = field
        self['level'] = level
        self['lead'] = lead
        self['months'] = (1,1)
        self['icdates'] = self._set_icdates()
        slead = "%03d" % self['lead']
        self.json = '_'.join(['spec', field, level, stream, centre, ftype, self['version'], slead]) + '.json'

    def compute(self):
        iclist = [int(datetime.strftime(ic, '%Y%m%d')) for ic in self['icdates']]
        self['en'] = []
        self['enrot'] = []
        self['endiv'] = []
        self['deg'] = []
        for m in range(min(self['months']), max(self['months'])+1):
            print("Running month "+str(m))
            month = "%02d" % m

            # Retrieve matching fields
            if self['field'] == 'ke':
                fld = self._read_fcst(month, iclist, fld='u')
                df = fld.data()[:,::-1,:]
                fld = self._read_fcst(month, iclist, fld='v')
                df2 = fld.data()[:,::-1,:]
            else:
                fld = self._read_fcst(month, iclist)
                df = fld.data()[:,::-1,:]
                df2 = None
                
            # Setup for grid-to-grid interpolation
            dims = df[0].shape
            lat0 = df[0,0,0]
            dlat = (df[0,dims[0]-1,0] - lat0) / (dims[0]-1)
            lon0 = df[1,0,0]
            dlon = (df[1,0,dims[1]-1] - lon0) / (dims[1]-1)
            sgrid = rmn.defGrid_L(dims[1], dims[0], lat0, lon0, dlat, dlon)
            dgrid = rmn.defGrid_G(dims[1]+1, dims[0])
            gridset = rmn.ezdefset(dgrid['id'], sgrid['id'])

            # Compute spherical harmonics for each grid
            for i in range(2, df.shape[0]):
                try:
                    (fout, fout2) = rmn.ezuvint(dgrid['id'], sgrid['id'],
                                                df[i,:,:].T, df2[i,:,:].T)
                    (deg, en, en_rot, en_div) = \
                        spharm.do_transform_ke(fout, fout2, assume_latlon='lonlat')
                    self['enrot'].append(en_rot.tolist())
                    self['endiv'].append(en_div.tolist())
                except TypeError:
                    fout = rmn.ezsint(dgrid['id'], sgrid['id'], df[i,:,:].T)
                    (deg, en) = \
                        spharm.do_transform_scalar(fout, assume_latlon='lonlat')
                if len(self['deg']) == 0: self['deg'] = deg.tolist()
                self['en'].append(en.tolist())
                            
    def write(self, fn):
        filt_dict = {
            key: value
            for key, value in self.items()
            if key not in ["an", "icdates"]
        }
        with open(fn, "w") as fd:
            json.dump(filt_dict, fd, indent=4)

    def _read_fcst(self, month, iclist, fld=None):
        rfld = fld or self['field']
        fname = '_'.join([rfld, self['level'], self['stream'], self['centre'],
                          self['ftype'], self['version'], 'pl', month])+'.grib2'
        fpath = os.path.join('forecasts', self['stream'], fname)
        ds = ek.data.from_source('file', fpath)
        flds = ek.data.SimpleFieldList()
        for fld in ds.sel(step=self['lead']):
            if fld.metadata("dataDate") in iclist:
                flds.append(fld)
        return(flds)

    def _set_icdates(self):
        datelist = [datetime.strptime('2024010100', '%Y%m%d%H')]
        while int(datetime.strftime(datelist[-1], '%Y')) == 2024:
            datelist.append(datelist[-1] + timedelta(days=3))
        return(datelist)

                
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('field', help='Field to score (ke for kinetic energy)')
    parser.add_argument('level', help='Pressure level to score (hPa)')
    parser.add_argument('-c', '--centre', default='cwao', help='WMO centre identifier (4-character)')
    parser.add_argument('-s', '--stream', default='oic', help='WP-MIP stream (oic/sic)')
    parser.add_argument('-t', '--type', default='pm', help='Model type (pm/hy/ai)')
    parser.add_argument('-r', '--version', default=0, help='Model version (integer)')
    parser.add_argument('-l', '--lead', default=120, help='Forecast lead time (h)')
    args = parser.parse_args()

    sp = Spec(args.centre, args.stream, args.type, args.version, args.field, args.level, args.lead)
    sp.compute()
    sp.write(sp.json)
