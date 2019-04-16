#!/usr/bin/env python3
# Provide a source copy of this setup along with a copy of the fargo3d source code
# and apply the necessary patches.
import argparse
import os
from subprocess import run, PIPE
import tempfile
import uuid

setupName = "hd163296"

def main():

    ### Command line arguments parser
    parser = argparse.ArgumentParser()
    parser.add_argument('--srcdir',
                        default=os.path.expanduser('~/repo/fargo3d/multifluid'),
                        help='Directory where the Fargo3d code is stored')
    parser.add_argument('--outdir', help='Directory to copy the code to')
    args = parser.parse_args()

    # directory where the setup and scripts are
    setupDir = os.path.dirname(argparse._sys.argv[0])

    outDir = provide(args.srcdir, setupDir, outDir=args.outdir)
    print(outDir)

def provide(srcDir, setupDir, outDir=None):
    if outDir is None: #provide a temp dir to work in
        outDir = tempfile.mkdtemp(prefix=setupName)
    else:
        outDir = os.path.abspath(outDir)

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
        copy( [os.path.join(srcDir, s) for s in [
            'Makefile',
            'scripts',
            'src',
            'license.txt',
            'std',
            'test_suite',
            'utils',
            'setups',
            'planets'
        ]], dst)

    copy_fargo3d(outDir)

    os.makedirs( os.path.join(outDir, 'setups'), exist_ok=True)
    os.makedirs( os.path.join(outDir, 'outputs'), exist_ok=True)

    ### copy the setup files
    #copy( os.path.join(setupDir, setupName),
    #      os.path.join(outDir, 'setups'))

    ### copy runscripts
    copy( os.path.join(setupDir, 'job'), outDir)

    ### copy patches
    copy( os.path.join(setupDir, 'patches'), outDir)

    ### link the run script
    os.symlink( 'job/run.sh', os.path.join(outDir, 'run.sh'))

    ### give the job a uuid
    with open(os.path.join(outDir, 'job/uuid.txt'), 'w') as idfile:
        idfile.write(str(uuid.uuid4()))

    return outDir

if __name__=="__main__":
    main()
