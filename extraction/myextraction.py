#!/usr/bin/env python

"""
Manual extraction
Modified to be able to run without issues on makeDirArrColumn() called on a AntennaPairStMan
It also allows to download selected fields (and not all the fields) that cover a cluster

G. Di Gennaro
Feb 2026
"""

from __future__ import print_function
import argparse
import sys
from run_extraction_pipeline import create_ds9_region,do_run_extract
from reprocessing_utils import prepare_field
import sys
from astropy.coordinates import SkyCoord,get_icrs_coordinates
import astropy.units as u
from auxcodes import run,report,warn,die,separator
import requests
import json
from astropy.table import Table
import glob
import os


# Call with the following signature: object name, [image size in deg], [ra,dec].
# If ra and dec are not provided the code will attempt to work them out from the name

parser = argparse.ArgumentParser(description='Run extraction of a cluster in LoTS')
parser.add_argument('-i','--clustername', help='cluster name, if you want to extract a single cluster', default='', required=False, type=str)
parser.add_argument('--RA', help='cluster RA (in deg)', required=False, type=float)
parser.add_argument('--DEC', help='cluster DEC (in deg)', required=False, type=float)
parser.add_argument('--size', help='size of box region (in deg)', default=0.4, required=True, type=float)
parser.add_argument('--fields', nargs='+', help='List of fields to download.')
args = vars(parser.parse_args())

#if len(sys.argv)==1:
#  die('Call this code with an object name and optional field size (in deg) and RA, DEC (in deg)',database=False)

separator('Finding target name and position')

target=args['clustername']
try:
  size=float(args['size'])
except:
  size=0.5

ra=None
dec=None
try:
  ra=float(args['RA'])
  dec=float(args['DEC'])
except:
    pass

print('Inputs to extraction. Size %s, RA %s, Dec %s, Name %s'%(size,ra,dec,target))

if ra is None:
  if 'ILTJ' in target:
    s=target[4:]
    coord=s[0:2]+':'+s[2:4]+':'+s[4:9]+' '+s[9:12]+':'+s[12:14]+':'+s[14:]
    sc = SkyCoord(coord,unit=(u.hourangle,u.deg))
    ra=sc.ra.value
    dec=sc.dec.value
    print('Parsed coordinates to ra=%f, dec=%f' % (ra,dec))
  else:
    sc=get_icrs_coordinates(target)
    ra=sc.ra.value
    dec=sc.dec.value
    print('Coordinate lookup gives ra=%f, dec=%f' % (ra,dec))
    target.replace(' ','')
else:
  sc=SkyCoord(ra*u.deg,dec*u.deg)
        
separator('Getting pointing positions')

r=requests.get('https://lofar-surveys.org/static/lotss_aladin/pointings_db.json')
d=json.loads(r.text)

if not args['fields']:
  names=[]
  ras=[]
  decs=[]
  for e in d:
    if e[3]=='Done':
      names.append(str(e[0]))
      ras.append(e[1])
      decs.append(e[2])

  t=Table(data=[names,ras,decs],names=['Field','ra','dec'])
  fsc=SkyCoord(t['ra']*u.deg,t['dec']*u.deg)
  t['sep']=sc.separation(fsc)

  fields=t[t['sep']<2.2*u.deg]
else:
  fields = args['fields']

if len(fields)==0:
  die('No fields within 2.2 degrees of pointing position',database=False)

print('We will use the following fields')
print(fields)

separator('Making working directory')

startdir = os.getcwd()
if not os.path.isdir(target):
  os.mkdir(target)
os.chdir(target)

if not os.path.exists('%s.reg'%(target)): # in this way, if we want to reprocess with a different box it does not create it
  create_ds9_region('%s.ds9.reg'%(target),ra,dec,size) # this needs to be modified for a more intelligent DS9 region based on flux  

separator('Downloading field data')

for f in fields:
  try:
    field=f['Field']
  except:
    field = f
  report('Doing field '+field)
  fdir=startdir+'/'+target+'/'+field
  if os.path.isdir(fdir):
    if len(glob.glob(fdir+'/*.ms.archive'))>0:
      warn('Field directory already contains MSs, skipping download')
      continue
  prepare_field(field,fdir,verbose=True)


separator('Running subtraction')

for f in fields:
  try:
    field=f['Field']
  except:
    field = f
  fdir=startdir+'/'+target+'/'+field
  os.chdir(fdir)
  
  executionstr = 'python '+os.environ['EXTRACTION_PATH']+'/sub-sources-outside-region.py -b ../%s.ds9.reg -p %s'%(target,target) 
  run(executionstr,database=False)

separator('Move subtracted datasets to working directory')

wd=startdir+'/'+target
run('cd %s; mv */*.dysco.sub.shift.avg.weights.ms.archive? .' %wd, database=False)

separator('Done!')

