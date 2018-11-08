#!/usr/bin/env python3
# Provide a source copy of this setup along with a copy of the fargo3d source code
# and apply the necessary patches.
import argparse
import os
from subprocess import run, PIPE
import tempfile
import uuid

setupName = "TD2P"

### Command line arguments parser
parser = argparse.ArgumentParser()
parser.add_argument('--fargo3d-src',
                    default=os.path.expanduser('~/repo/fargo3d/fargo3d-1.3'),
                    help='Directory where the Fargo3d code is stored')
parser.add_argument('--name', help='Name of the project/job')
parser.add_argument('--outdir', help='Directory to copy the code to')
args = parser.parse_args()

fargoDir = args.fargo3d_src
if args.outdir is None: #provide a temp dir to work in
    outDir = tempfile.mkdtemp(prefix=setupName)    
else:
    outDir = os.path.abspath(args.outdir)
# place the files into its own folder inside the outdir if name is given
if args.name is not None:
    outDir = os.path.join( outDir, args.name)
    os.makedirs(outDir)
print(outDir)
# directory where the setup and scripts are
setupDir = os.path.dirname(argparse._sys.argv[0])

### Make sure the output dir exists
if not os.path.exists(outDir):
    os.makedirs(outDir)

### copy the necessary fargo3d sources
def copy(src, dst):
    if isinstance(src, str):
        src = [src]
    run(['cp', '-r'] + src + [dst])

def copy_fargo3d(dst):
    # copy the fargo3d src to the output dir
    copy( [os.path.join(fargoDir, s) for s in [
        'Makefile',
        'scripts',
        'src',
        'license.txt',
        'std',
        'test_suite',
        'utils'
    ]], dst)

copy_fargo3d(outDir)
    
os.makedirs( os.path.join(outDir, 'setups'), exist_ok=True)
os.makedirs( os.path.join(outDir, 'outputs'), exist_ok=True)

### copy the setup files
copy( os.path.join(setupDir, setupName),
      os.path.join(outDir, 'setups'))

### copy runscripts
copy( os.path.join(setupDir, 'job'), outDir)

### copy patches
copy( os.path.join(setupDir, 'patches'), outDir)

### link the run script
os.symlink( 'job/run.sh', os.path.join(outDir, 'run.sh'))

### give the job a uuid
with open(os.path.join(outDir, 'job/uuid.txt'), 'w') as idfile:
    idfile.write(str(uuid.uuid4()))
