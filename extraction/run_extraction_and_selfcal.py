import argparse
import glob, os, sys
import numpy as np
import pickle
try:
    import bdsf as bdsm
except ImportError:
    import lofar.bdsm as bdsm
from crossmatch_utils import *
import pyregion

import time
start_time = time.time()

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy import wcs
from astropy.wcs import WCS
from auxcodes import separator


"""
modify these paths to link to the macaroon and to the location of the extraction folder
"""
os.environ['RCLONE_CONFIG_DIR'] = './'
os.environ['EXTRACTION_PATH'] = './'


from reprocessing_utils import *

DATADIR = './'

def upload_extract(cname,uploadfilename):
	rclonepath = os.environ['RCLONE_CONFIG_DIR']
	print('rclone  --multi-thread-streams 1 --config=%s/maca_sksp_disk_extract.conf copy %s maca_sksp_disk_extract:%s/selfcal/'%(rclonepath,uploadfilename,cname))
	os.system('rclone  --multi-thread-streams 1 --config=%s/maca_sksp_disk_extract.conf copy %s maca_sksp_disk_extract:%s/selfcal/'%(rclonepath,uploadfilename,cname))


parser = argparse.ArgumentParser(description='Run extraction and selfcalibration of clusters in LoTSS; you can give either a catalog (FITS format) or the cluster name')
parser.add_argument('--doextraction', help='set it True if you want to extract clusters from archive', action='store_true')
parser.add_argument('--doselfcal', help='set it True if you want to run facetselfcal', action='store_true')
parser.add_argument('-i','--clustername', help='cluster name, if you want to extract a single cluster', default='', required=False, type=str)
parser.add_argument('--RA', help='cluster RA (in deg)', required=False, type=float)
parser.add_argument('--DEC', help='cluster DEC (in deg)', required=False, type=float)
parser.add_argument('--size', help='size of box region (in deg)', default=0.4, required=False, type=float)
parser.add_argument('--fields', nargs='+', help='For a single target, list of fields to download', required=False)
parser.add_argument('-c','--catalog', help='Catalog to use from which extract clusters', required=False, type=str)


args = vars(parser.parse_args())

if args['catalog'] and args['clustername']:
  print ("Error: either give a single target name or a cluster catalog")
  sys.exit()

elif args['catalog'] and not args['clustername']:
  print ("use catalog:", args['catalog'])
  data = fits.open(args['catalog'])[1].data
  clusterlist = np.array(data['Name']) 
  ra         = np.array(data['RAJ2000'])
  dec        = np.array(data['DEJ2000'])

elif args['clustername'] and not args['catalog']:
  clusterlist = [args['clustername']]
  try:
    ra         = [args['RA']]
    dec        = [args['DEC']]
  except:
    ra, dec = [None], [None]
else:
  print ("Error: give a single target name or a cluster catalog.")
  sys.exit()


for i, cluster in enumerate(clusterlist):
  print (os.getcwd())
  name  = cluster.replace(' ','')
  RA    = ra[i]
  DEC   = dec[i]
  print ("CLUSTER:", name, RA, DEC)
  
  if not args['size']: 
    size = 0.4
  else:
    size = args['size']

  if args['doextraction']:   
    if len(glob.glob(DATADIR+name+"/"+"P???+??*.dysco.sub.shift.avg.weights.ms.archive*")) == 0:
      #cmd = 'python '+os.environ['EXTRACTION_PATH']+'/myextraction.py %s %s %s %s'%(name,size,RA,DEC)
      try: #if len(args['fields']) > 0:
        cmd = 'python '+os.environ['EXTRACTION_PATH']+'/myextraction.py -i %s --size %s --RA %s --DEC %s --fields %s'%(name,size,RA,DEC," ".join(map(str,args['fields'])))
      except: #else:
         cmd = 'python '+os.environ['EXTRACTION_PATH']+'/myextraction.py -i %s --size %s --RA %s --DEC %s'%(name,size,RA,DEC)
      print (cmd)
      os.system(cmd)
			os.system('tar -zcf extractedMS.tar.gz *.dysco.sub.shift.avg.weights.ms.archive?') #save the original source-subtracted MSes
    else:
      print ("Extraction already done")
     
  else:
    print ("No extraction requested for %s"%name)
    if len(glob.glob(DATADIR+name+"/"+"P???+??*.dysco.sub.shift.avg.weights.ms.archive*")) == 0:
      os.system("mv "+DATADIR+name+"/*/"+"P???+??*.dysco.sub.shift.avg.weights.ms.archive* "+DATADIR+name+"/.")
  
  if os.path.exists(DATADIR+name) and args['doselfcal'] and len(glob.glob(DATADIR+name+"/"+"P???+??*.dysco.sub.shift.avg.weights.ms.archive*")) > 0:
    if not os.path.exists(DATADIR+name+'/'+name+'.ds9.tar.gz'):      
      # RUN SELF CALIBRATION AND COPYING IN THE SURFSARA DATABASE
      CLUSTERDIR = DATADIR+name+'/'
      os.chdir(CLUSTERDIR)
      print (os.getcwd())
      
      cmd  = 'python '+os.environ['EXTRACTION_PATH']+'/lofar_facet_selfcal/facetselfcal.py '
      cmd += '-b %s.ds9.reg '%name
      cmd += '--auto '  
			cmd += '--metadata-compression False '
      cmd += '-i %s *dysco.sub.shift.avg.weights.ms.archive?'%name
      
      print (cmd)
      os.system(cmd)
    
      os.chdir('../')

    else:
      print ("Selfacl already done")   

  else:
    print ("check extraction output") 
    
  timerun = round((time.time() - start_time)/3600,1)
  separator("%s hr"%(str(timerun)))



