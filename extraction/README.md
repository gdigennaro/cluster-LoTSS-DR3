The script `run_extraction_and_selfcal.py` works both on single clusters, with coordinates cluster name provided, and with a FITS table with a list of clusters; in this latter case, the table needs to have Name, RAJ2000, DEJ2000 columns). This is made of two independent steps:

- the extraction, where the pointing is downloaded and sources outside a given region are subtracted. It requires a macaroon to connect to the LoTSS database (provided only on request by gabriella.digennaro@inaf.it)
- the selfcal, performed using `facetselfcal v17.14` in auto mode (4 cycles of phase only and 6 cycles of amp+phase calibration). Possibly, you can add/decrease selfcal cycle by adding `--stop NCYCLES` to the `facetselfcal.py` command (see https://github.com/rvweeren/lofar_facet_selfcal/tree/main). I suggest to download always the most updated version of `facetselfcal.py` and to place it in `os.environ['EXTRACTION_PATH']` 
 
It needs to be ran inside `flocs v6.0` (https://public.spider.surfsara.nl/project/lofarvwf/fsweijen/containers/).

It requires the definition of `os.environ['RCLONE_CONFIG_DIR']` and `os.environ['EXTRACTION_PATH']` to the location of the macaroon and to the `extraction` folder (where `run_extraction_and_selfcal.py`, `myextraction.py` and `sub=sources-outside-region.py` are located). For single targets, one can also choose to download selected fields (parameter `[--fields]`)


From the run directory, it will create subfolder(s) with the CLUSTERNAME(s).

Examples of how to run the script, for a cluster list:

`python run_extraction_and_selfcal.py -c CLUSTERLIST.fits [–-doextraction] [--doselfcal]`

and for a single cluster:

`python run_extraction_and_selfcal.py -i CLUSTERNAME --RA CLUSTER_RA --DEC CLUSTER_DEC [--size 0.4] [–-doextraction] [--doselfcal] [--fields]`

Issues and fix are available at this doc: https://docs.google.com/document/d/1vPrjBnpE0didmAeQlA2b2FeRBY8AqJ5rVkOfs9CEp-g/edit?usp=sharing
