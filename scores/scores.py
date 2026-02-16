#!/usr/bin/env python
import glob, os, sys, json
import earthkit as ek
import xarray as xr
import pandas as pd
import numpy as np
from datetime import datetime, timedelta

class Score(dict):

    def __init__(self, centre, stream, ftype, version, field, level, analyses):
        self['centre'] = centre
        self['stream'] = stream
        self['ftype'] = ftype
        self['version'] = "%02d" % int(version)
        self['field'] = field
        self['level'] = level
        self['analyses'] = analyses
        self['months'] = (1,12)
        self['icdates'] = self._set_icdates()
        self['an'] = self._read_analyses(analyses)
        self.json = '_'.join([field, level, stream, centre, ftype, self['version']]) + '.json'

    def compute(self, score):
        compute_score = getattr(self, '_'+score)
        self[score] = {}
        for m in range(min(self['months']), max(self['months'])+1):
            print("Running month "+str(m))
            month = "%02d" % m
            fname = '_'.join([self['field'], self['level'], self['stream'], self['centre'],
                              self['ftype'], self['version'], 'pl', month])+'.grib2'
            fpath = os.path.join('forecasts', self['stream'], fname)
            ds = ek.data.from_source('file', fpath).to_xarray()
            for init, this_init in ds.groupby("forecast_reference_time"):
                init_dt = datetime.fromisoformat(str(this_init.forecast_reference_time.values[0]))
                if init_dt not in self['icdates']:
                    continue
                for fcst, this_fcst in this_init.groupby("step"):
                    hr = round(float(fcst / 3600e9))
                    if hr % 12 != 0 or hr > 240: continue
                    vtime = init + fcst
                    for an in self['an'].values():
                        fld_an = an.sel(forecast_reference_time=vtime)
                        scorev = compute_score(this_fcst, fld_an)
                    try:
                        self[score][hr].append(scorev)
                    except(KeyError):
                        self[score][hr] = [scorev]

    def write(self, fn):
        filt_dict = {
            key: value
            for key, value in self.items()
            if key not in ["an", "icdates"]
        }
        with open(fn, "w") as fd:
            json.dump(filt_dict, fd, indent=4)

    def _rmse(self, dsf, dsa):
        msf = xr.ufuncs.cos(xr.ufuncs.radians(dsa.latitude))
        rmse = (msf * (dsf - dsa)**2).mean(dim=['latitude', 'longitude']) / msf.mean()
        return(rmse[self['field']].values.item())

    def _bias(self, dsf, dsa):
        msf = xr.ufuncs.cos(xr.ufuncs.radians(dsa.latitude))
        bias = (msf * (dsf - dsa)).mean(dim=['latitude', 'longitude']) / msf.mean()
        return(bias[self['field']].values.item())
        
    def _read_analyses(self, analyses):
        anvals = {}
        an_data = []
        for analysis in self['analyses']:
            for m in range(min(self['months']), max(self['months'])+2):
                month = "%02d" % m
                fname = '_'.join([self['field'], self['level'], analysis, 'pl', month])+'.grib2'
                fpath = os.path.join('analyses', fname)
                ds = ek.data.from_source('file', fpath).to_xarray()
                an_data.append(ds)
            print("done reading")
            anvals[analysis] = (xr.concat(an_data, 'forecast_reference_time', data_vars="minimal",
                                          coords="minimal", compat="override"))
            print("done concat")
        return(anvals)

    def _set_icdates(self):
        datelist = [datetime.strptime('2024010100', '%Y%m%d%H')]
        while int(datetime.strftime(datelist[-1], '%Y')) == 2024:
            datelist.append(datelist[-1] + timedelta(days=3))
        return(datelist)

                
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('field', help='Field to score')
    parser.add_argument('level', help='Pressure level to score (hPa)')
    parser.add_argument('-c', '--centre', default='cwao', help='WMO centre identifier (4-character)')
    parser.add_argument('-s', '--stream', default='oic', help='WP-MIP stream (oic/sic)')
    parser.add_argument('-t', '--type', default='pm', help='Model type (pm/hy/ai)')
    parser.add_argument('-r', '--version', default=0, help='Model version (integer)')
    parser.add_argument('--score', action='append', help='Scores to compute')
    parser.add_argument('--anal', default=['ecmf'], action='append', help='Analyses for scoring')
    args = parser.parse_args()

    sc = Score(args.centre, args.stream, args.type, args.version, args.field, args.level, args.anal)
    for score in args.score:
        sc.compute(score)
    sc.write(sc.json)
        
