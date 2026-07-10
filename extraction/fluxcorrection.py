"""
This script aims to correct the flux scale to the LoTSS-DR3 mosaics.
"""


#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import urllib.request
import urllib.parse
import os
import re
import glob
import sys
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from regions import Regions

import warnings
warnings.filterwarnings("ignore")

# Try importing bdsf
try:
  import bdsf
except ImportError:
  print("Warning: 'bdsf' module not found. Please ensure PyBDSF is installed.")

def download_lofar_mosaic(cluster_name, ra, dec, size_arcmin, out_file='dr3-cutout.fits'):
  """Downloads a FITS cutout from the LOFAR DR3 server using the cluster name."""
  print(f"\n--- Downloading LOFAR DR3 Mosaic for '{cluster_name}' ---")
  
  try:
    encoded_pos = urllib.parse.quote(cluster_name)
    url = f"https://lofar-surveys.org/dr3-cutout.fits?pos={encoded_pos}&size={size_arcmin}"  
    print(f"Requesting: {url}")
    try:
      urllib.request.urlretrieve(url, out_file)
      print(f"Successfully downloaded to: {out_file}")
      return out_file
    except:
      encoded_pos = urllib.parse.quote(ra+' '+dec)
      url = f"https://lofar-surveys.org/dr3-cutout.fits?pos={encoded_pos}&size={size_arcmin}"
      urllib.request.urlretrieve(url, out_file)
      print(f"Successfully downloaded to: {out_file}")

  except Exception as e:
    print(f"Error downloading mosaic: {e}")
    raise

def run_pybdsf(fitsname, detectimage=None, out_reg='sources.reg', trim_frac=0.15):
  """Runs PyBDSF to detect sources and generate a raw DS9 region file."""
  print(f"\n--- Running PyBDSF Source Detection on {fitsname} ---")
  
  bdsf_params = {
      'thresh_isl': 4.0, 
      'thresh_pix': 5.0, 
      'rms_box': (150, 15), 
      'rms_map': True, 
      'mean_map': 'zero',
      'ini_method': 'intensity', 
      'adaptive_rms_box': True, 
      'adaptive_thresh': 150, 
      'rms_box_bright': (60, 15),
      'group_by_isl': False, 
      'group_tol': 10.0, 
      'output_opts': True, 
      'output_all': True, 
      'atrous_do': False,      # Essential: prevents heavy fitting of extended emission
      'flagging_opts': True, 
      'flag_maxsize_fwhm': 0.5,
      'advanced_opts': True, 
      'blank_limit': None
  }
  
  # Calculate and apply the trim_box if requested
  if trim_frac > 0:
    with fits.open(fitsname) as hdul:
      ny, nx = hdul[0].data.shape[-2:]
        
    xmin = int(nx * trim_frac)
    xmax = int(nx * (1.0 - trim_frac))
    ymin = int(ny * trim_frac)
    ymax = int(ny * (1.0 - trim_frac))
    
    print(f"Applying trim_box to avoid outer {trim_frac*100}% of image edges.")
    print(f"Searching only within: X({xmin}:{xmax}), Y({ymin}:{ymax})")
    bdsf_params['trim_box'] = (xmin, xmax, ymin, ymax)
  
  if detectimage is not None:
    bdsf_params['detection_image'] = detectimage
      
  img = bdsf.process_image(fitsname, **bdsf_params)
  img.write_catalog(outfile=out_reg, catalog_type='srl', format='ds9', clobber=True)
  return out_reg

def get_beam_area_pixels(header):
  """Calculates the beam area in pixels from FITS header safely."""
  bmaj = header.get('BMAJ', 0) * 3600
  bmin = header.get('BMIN', 0) * 3600
  
  if 'CDELT1' in header:
    pixscale = abs(header['CDELT1'] * 3600.)
  elif 'CD1_1' in header:
    pixscale = abs(header['CD1_1'] * 3600.)
  else:
    pixscale = 1.0
      
  if bmaj == 0 or bmin == 0: return 1.0
  return (np.pi * bmaj * bmin) / (4 * np.log(2) * pixscale**2)

def process_image(img_path, source_reg_path, flux_scale_err=0.10):
  """Extracts fluxes and stores them in a dictionary keyed by region index."""
  print(f"Extracting fluxes from: {img_path}")
  with fits.open(img_path) as hdul:
    data = hdul[0].data.squeeze()
    header = hdul[0].header
    wcs = WCS(header).celestial
      
  beam_area_pix = get_beam_area_pixels(header)
  stats = {'sources': {}} 

  source_regions = Regions.read(source_reg_path, format='ds9')
  
  for i, reg in enumerate(source_regions):
    pix_reg = reg.to_pixel(wcs)
    reg_mask = pix_reg.to_mask()
    
    if reg_mask is None: continue
    img_mask = reg_mask.to_image(data.shape)
    if img_mask is None: continue
    
    source_pixels = data[img_mask.astype(bool)]
    source_pixels = source_pixels[~np.isnan(source_pixels)]
    
    if len(source_pixels) == 0: continue

    flux_jy = np.sum(source_pixels) / beam_area_pix
    peak_jy = np.max(source_pixels)
    
    src_stats = {
        'flux_jy': flux_jy,
        'peak_jy': peak_jy,
        'total_flux_error_jy': flux_scale_err * abs(flux_jy)
    }
    stats['sources'][i] = src_stats

  return stats

def create_comparison_plot(df, cluster_name, label1="Self-Cal Image", label2="DR3 Mosaic", out_png='flux_comparison.png'):
  """Generates a log-scale 1:1 comparison scatter plot."""
  print("\n--- Generating Comparison Plot ---")
  
  df['Flux_Diff'] = df['Flux_img1'] - df['Flux_img2']
  df['Flux_Ratio'] = np.where(df['Flux_img2'] != 0, df['Flux_img1'] / df['Flux_img2'], np.nan)

  print(f"Mean Difference ({label1} - {label2}): {df['Flux_Diff'].mean():.6e} Jy")
  print(f"Median Difference: {df['Flux_Diff'].median():.6e} Jy")
  print(f"Mean Ratio ({label1} / {label2}): {df['Flux_Ratio'].mean():.4f}")
  print(f"Median Ratio: {df['Flux_Ratio'].median():.4f}")

  plt.figure(figsize=(6, 6))
  
  plt.errorbar(df['Flux_img2'], df['Flux_img1'], 
                xerr=np.abs(df['Flux_err_img2']), yerr=np.abs(df['Flux_err_img1']), 
                fmt='o', alpha=1, color='k', ecolor='gray', capsize=0, label='extracted point sources')
      
  min_val = min(df['Flux_img2'].min(), df['Flux_img1'].min())
  max_val = max(df['Flux_img2'].max(), df['Flux_img1'].max())
  min_plot = min_val * 0.8 if min_val > 0 else min_val - 0.01
  max_plot = max_val * 1.2

  offset = df['Flux_Ratio'].median()
  
  plt.plot([min_plot, max_plot], [min_plot, max_plot], color='k', linestyle='--', linewidth=1, label='1:1 ratio (perfect match)')
  plt.plot([min_plot, max_plot], [offset*min_plot, offset*max_plot], color='cornflowerblue', linestyle='-.', linewidth=2, label='corrected ratio (median=%s)'%(round(offset,2)))
  
  if df['Flux_img2'].min() > 0 and df['Flux_img1'].min() > 0:
    plt.xscale('log')
    plt.yscale('log')
    plt.title(f'{cluster_name} - Point Source Flux Comparison', fontsize=12)
  else:
    plt.title(f'{cluster_name} - Point Source Flux Comparison', fontsize=12)

  plt.xlabel(f'{label2} Flux (Jy)', fontsize=11)
  plt.ylabel(f'{label1} Flux (Jy)', fontsize=11)
  plt.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, right=True, left=True, labeltop=False, labelright=False)
  #plt.grid(True, which="both", ls="--", alpha=0.4)
  l = plt.legend(fontsize=8)
  plt.setp(l.texts)

  
  if not (df['Flux_img2'].min() > 0 and df['Flux_img1'].min() > 0): plt.axis('equal')
      
  plt.xlim(min_plot, max_plot)
  plt.ylim(min_plot, max_plot)
  #plt.tight_layout()

  plt.savefig(out_png, dpi=300, bbox_inches='tight')
  print(f"Plot saved successfully to: {out_png}")


def remove_text_from_region_file(input_path, output_path):
  if output_path is None:
    output_path = input_path

  # Matches  text={  ...  }  where the value may contain nested braces.
  # The pattern uses a simple repeated-group approach that handles one level
  # of nesting (sufficient for DS9 region labels).
  text_pattern = re.compile(r'\s*text=\{[^{}]*(?:\{[^{}]*\}[^{}]*)?\}')

  with open(input_path, "r") as fh:
    lines = fh.readlines()

  cleaned_lines = []
  for line in lines:
    if line.lstrip().startswith("global"):
      # Remove text={...} token(s)
      line = text_pattern.sub("", line)
      # Collapse multiple consecutive spaces (but preserve leading spaces)
      line = re.sub(r'(?<=\S) {2,}', ' ', line)
    cleaned_lines.append(line)

  cleaned_content = "".join(cleaned_lines)

  with open(output_path, "w") as fh:
    fh.write(cleaned_content)

  print(f"Cleaned region file written to: {output_path}")
  return cleaned_content


parser = argparse.ArgumentParser(description="Pipeline Sanity Check: DR3 Download, PyBDSF Detection, Point Source Filtering & Comparison")
parser.add_argument('--catalog', help='Catalog to use from which extract clusters', required=False, type=str)  
parser.add_argument("-i1", "--image1", required=False, help="Path to primary (reference) FITS image")
parser.add_argument("-i2", "--image2", required=False, help="Path to secondary (either your image or to download from DR3) FITS image")
parser.add_argument("-ra", "--ra", required=False, help="ra of the cluster (in case it is not recognised by name")
parser.add_argument("-dec", "--dec", required=False, help="dec of the cluster (in case it is not recognised by name")
parser.add_argument("-c", "--clustername", required=False, help="Name of the cluster (e.g., 'Abell 2256')")
parser.add_argument("-s", "--size", type=float, default=60.0, help="Size of the LOFAR cutout to download in arcminutes (default: 60)")
parser.add_argument("-d", "--detect_image", default=None, help="Path to detection image for PyBDSF (Optional)")
parser.add_argument("-fe", "--flux_err", type=float, default=0.10, help="Fractional flux scale error (default: 0.10)")
parser.add_argument("-pt", "--point_thresh", type=float, default=1.5, help="Max Integrated/Peak flux ratio to be a point source (default: 1.5)")
parser.add_argument("-t", "--trim", type=float, default=0.15, help="Fraction of the image edges to trim to avoid noise (default: 0.15 for 15%)")
parser.add_argument("-o", "--out_prefix", default="sanity_check", help="Optional extra prefix for output files")
parser.add_argument("--start", help="if running pybdsf", action='store_true')
parser.add_argument("--docleanup", help="if running pybdsf", action='store_true')
args = parser.parse_args()


if args.catalog and args.clustername:
  print ("Error: either give a single target name or a cluster catalog")
  sys.exit()

elif args.catalog and not args.clustername:
  print ("use catalog:", args.catalog)
  data = fits.open(args.catalog)[1].data
  clusterlist = np.array(data['Name']) 
  ra         = np.array(data['RAJ2000'])
  dec        = np.array(data['DEJ2000'])

elif args.clustername and not args.catalog:
  clusterlist = [args.clustername]
  try:
    ra         = [args.ra]
    dec        = [args.dec]
  except:
    ra, dec = [None], [None]
else:
  print ("Error: give a single target name or a cluster catalog.")
  sys.exit()

fluxcorr = np.array([])
for i, cluster in enumerate(clusterlist):
  # Create a clean, filesystem-safe version of the cluster name
  safe_cluster = cluster.replace(" ", "").replace(",", "")

  if args.image1: 
    image1 = args.image
  else:
    image1 = '../%s/LOFAR/%s_009-MFS-image.fits'%(safe_cluster,safe_cluster)
  
  # Construct a base name that strictly includes the cluster name
  base_name = f"{safe_cluster}_{args.out_prefix}"

  if not args.image2:
    # Download DR3 Mosaic 
    dr3_mosaic = f"{base_name}_DR3-cutout.fits"
    download_lofar_mosaic(safe_cluster, ra[i], dec[i], args.size, out_file=dr3_mosaic)
  else:
    dr3_mosaic = args.image2
  
  # Run PyBDSF to get raw regions
  raw_region_file = f"{base_name}_raw_sources.reg"
  if args.start:
    run_pybdsf(image1, detectimage=args.detect_image, out_reg=raw_region_file, trim_frac=args.trim)
  
  # Extract fluxes
  print("\n--- Calculating Fluxes ---")
  stats1 = process_image(image1, raw_region_file, args.flux_err)
  stats2 = process_image(dr3_mosaic, raw_region_file, args.flux_err)

  # Filter for Point Sources & Build Data
  print(f"\n--- Filtering Extended Sources (Int/Peak Ratio < {args.point_thresh}) ---")
  source_regions = Regions.read(raw_region_file, format='ds9')
  valid_regions = []
  data = []
  
  for i in range(len(source_regions)):
    # Make sure the source was successfully extracted in both images
    if i in stats1['sources'] and i in stats2['sources']:
      s1 = stats1['sources'][i]
      s2 = stats2['sources'][i]
      
      flux_img1 = s1['flux_jy']
      peak_img1 = s1['peak_jy']
      
      # Discard negative fluxes
      if flux_img1 <= 0 or peak_img1 <= 0: continue
      
      int_peak_ratio = flux_img1 / peak_img1
      
      # POINT SOURCE CHECK: Ratio must be between 0.5 and the threshold (1.5)
      if 0.8 <= int_peak_ratio <= args.point_thresh:
        valid_regions.append(source_regions[i])
        data.append({
            'ID': i + 1,
            'Int_Peak_Ratio': round(int_peak_ratio, 3),
            'Flux_img1': flux_img1,
            'Flux_err_img1': s1['total_flux_error_jy'],
            'Flux_img2': s2['flux_jy'],
            'Flux_err_img2': s2['total_flux_error_jy']
        })
              
  if len(source_regions) == 0:
    print("Warning: PyBDSF found 0 sources. Consider checking your threshold parameters or image quality.")
  else:
    print(f"Discarded {len(source_regions) - len(valid_regions)} extended sources/artifacts.")
    print(f"Kept {len(valid_regions)} true point sources.")
  
  # Save Filtered Region File for DS9 Validation
  filtered_reg_file = f"{base_name}_clean_point_sources.reg"
  with open(filtered_reg_file, 'w') as f:
    f.write("global color=green dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n")
    for reg in valid_regions:
      f.write(reg.serialize(format='ds9') + '\n')
  print(f"Saved clean point-source region file to: {filtered_reg_file}")
  remove_text_from_region_file(filtered_reg_file, filtered_reg_file)
  
  # Save CSV & Plot (only if we found sources)
  if len(valid_regions) > 0:
    df = pd.DataFrame(data)
    csv_file = f"{base_name}_fluxes.csv"
    df.to_csv(csv_file, index=False)
    print(f"Flux catalog saved to: {csv_file}")

    plot_file = f"{base_name}_plot.pdf"
    create_comparison_plot(df, safe_cluster, label1="self-cal", label2="DR3", out_png=plot_file)
  else:
    print("Not enough point sources survived the filtering to create a plot.")
  
  
  if docleanup:
    print("\n--- Cleaning Up Temporary Files ---")
    if os.path.exists(dr3_mosaic):
      os.remove(dr3_mosaic)
      print(f"Deleted downloaded DR3 cutout: {dr3_mosaic}")
    if os.path.exists(raw_region_file):
      os.remove(raw_region_file)
      print(f"Deleted raw un-filtered regions: {raw_region_file}")
  
  fluxcorr = np.append(fluxcorr, 1./df['Flux_Ratio'].median())
  print("\nSanity Check Pipeline Completed!")


if args.catalog and not args.clustername:
  updatedtab = args.catalog
  table = Table.read(updatedtab)
  table['fluxcorr'] = fluxcorr
  table.write(updatedtab, overwrite=True)
