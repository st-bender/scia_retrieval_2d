#!/usr/bin/env python
# vim: set fileencoding=utf-8

from __future__ import absolute_import, division, print_function

import os
import sys
import numpy as np
import datetime as dt
import dateutil.parser as dp
import scipy.interpolate as sip
import bisect
import json
try:
	from StringIO import StringIO
except ImportError:
	from io import BytesIO as StringIO

mkdate = lambda t: dt.datetime.strptime(bytes(t).decode(), "%Y-%m-%d")
mktime = lambda t: dt.datetime.strptime(bytes(t).decode(), "%H:%M")

def get_solar_data(fname, scale=1.0):
	sol_name = os.path.basename(fname).split('_')[1]
	sol_data = np.genfromtxt(fname, dtype=None,
			converters = {0: mkdate, 1: mktime})
	sol_dates = sol_data['f0']
	sol_values = sol_data['f2']*scale
	return sol_name, sol_dates, sol_values

def time_fun_vec(t, N=[1.0, 2.0], solvalues=[], offset=True, linear=False):
	w = 2.*np.pi/365.25
	try:
		sval = [s[t] for s in solvalues]
	except:
		sval = [np.median(np.asarray(s.values())) for s in solvalues]
	offs = [1.0] if offset else []
	lin = [t] if linear else []
	ret = offs + lin + [f(n*w*t) for f in [np.cos, np.sin] for n in N] + sval
	return ret

def time_fun_sol(tday, t0, solnames, solscales, sol_lags, solsubmin=1):
	sdicts = []
	for solname, solscale, sollag in zip(solnames, solscales, sol_lags):
		sname, sdates, svalues = get_solar_data(solname, solscale)
		if solsubmin > 0:
			if solsubmin == 1:
				svalues -= svalues.min()
			else:
				dds = [days, ds][solsubmin - 2]
				sidx = [list(sdates).index(t0 + dt.timedelta(x))
						for x in [(d-t0).days for d in sdates if (d-t0).days in dds]]
				svalues -= svalues[sidx].min()
		sdicts.append(dict(zip([(t - t0).days + sollag for t in sdates], svalues)))
	return time_fun_vec(tday, solvalues=sdicts)

def coeffs_alt_lat(coeffs, fields, alt=100.0, lat=67.5):
	lats = np.unique(coeffs['lats'])
	lats_right = lats + np.mean(0.5*np.diff(lats))
	j = bisect.bisect_left(lats_right, lat)
	if j >= lats.shape[0]:
		# decrease index if it is too large
		j = lats.shape[0] - 1
	l1 = lats[j]
	cfl1 = coeffs[coeffs['lats']==l1]
	alts = cfl1['alts']
	return (alts, cfl1[fields].copy().view('<f8').reshape((-1, len(fields))))

def main():
	ff = open(sys.argv[1], 'r')
	NO_json = json.load(ff)
	files = NO_json["files"]
	config = NO_json["config"]

	model_t0 = dt.datetime(2003, 7, 1)

	coeff_fname = files["data"]["regress"]
	data_name = config["dataname"].encode("utf-8")
	date = dp.parse(sys.argv[2])
	t_day = (date - model_t0).days

	try:
		alt = np.genfromtxt(StringIO(sys.argv[3].encode()), dtype='<f8', delimiter=",")
		lat = float(sys.argv[4])
	except:
		alt = 100.0
		lat = 67.5

	#print >> sys.stderr, alt, lat

	#basepath = os.path.dirname(sys.argv[0])
	#if basepath == '':
	#	basepath = '.'

	sol_files = files["solar"]
	solnames = config.get("solnames", ["lya", "kp"])
	#sol_data = [basepath + '/' + sol_files[sn] for sn in solnames]
	sol_data = [sol_files[sn] for sn in solnames]
	solscales = config.get("solscales", [1.0, 1.0])
	sol_lags = config.get("sol_lags", [0, 1])
	solsubmin = config.get("solsubmin", 1)

	fields = ["offset", "cos1", "sin1", "cos2", "sin2"] + solnames
	data = np.genfromtxt(coeff_fname, dtype=None, names=True)
	coeffs = data[data['name']==data_name]

	#print >> sys.stderr, date, coeff_fname, t_day
	#print >> sys.stderr, sol_files, solnames, sol_data

	tfv = np.asarray(time_fun_sol(t_day, model_t0, sol_data, solscales,
			sol_lags, solsubmin))
	alts, cvs = coeffs_alt_lat(coeffs, fields, lat=lat)
	mod_val = np.dot(cvs, tfv)
	mod_val_f = sip.InterpolatedUnivariateSpline(alts, mod_val, k=1)
	#mod_val_f = sip.UnivariateSpline(alts, mod_val, k=1)
	#mod_val_f = sip.interp1d(alts, mod_val, bounds_error=False, fill_value=0)
	print(' '.join(map(str, mod_val_f(np.atleast_1d(alt)))))

if __name__ == "__main__":
	main()
