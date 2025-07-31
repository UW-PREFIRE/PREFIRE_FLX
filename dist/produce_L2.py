"""
PROC_MODE = 1: Create part of a 2B-FLX product granule (using only AUX-MET)
PROC_MODE = 2: Create part of a 2B-FLX product granule (using AUX-MET + AUX-SAT)

This program requires python version 3.6 or later, and is importable as a 
python module.
"""

  # From the Python standard library:
from pathlib import Path
import os
import sys
import argparse
import subprocess

  # From other external Python packages:

  # Custom utilities:


#--------------------------------------------------------------------------
def main(use_saved_exe, interactive_MatLab):
    """Driver routine."""

    package_top_Path = Path(os.environ["PACKAGE_TOP_DIR"])
    
    sys.path.append(str(package_top_Path / "source" / "python"))

    this_environ = os.environ.copy()

    proc_mode = int(this_environ["PROC_MODE"])

      # Default product_fullver:
    if "PRODUCT_FULLVER" not in this_environ:
        this_environ["PRODUCT_FULLVER"] = "R01_P00"
    elif len(this_environ["PRODUCT_FULLVER"].strip()) == 0:
        this_environ["PRODUCT_FULLVER"] = "R01_P00"

    # Requires MatLab

    MatLab_src_dir = str(package_top_Path / "source" / "matlab")
    app_MatLab_prefix = "app_matlab"

    dist_Path = package_top_Path / "dist"

    if interactive_MatLab:
        this_environ["MATLABPATH"] = f"{MatLab_src_dir}:{dist_Path}"
        cmd = ["matlab", "-singleCompThread", "-nodisplay"]
    else:
        if use_saved_exe:
            cmd = [str(dist_Path / f"{app_MatLab_prefix}_run")]
        else:
            this_environ["MATLABPATH"] = f"{MatLab_src_dir}:{dist_Path}"
            cmd = ["matlab", "-singleCompThread", "-nodisplay", "-batch",
                   f"{app_MatLab_prefix}_run"]

    subprocess.run(cmd, env=this_environ)


if __name__ == "__main__":
    # Process arguments:
    arg_description = ("PROC_MODE = 1: Create part of a 2B-FLX product "
                       " granule (using AUX-MET surface type only)\n"
                       "PROC_MODE = 2: Create part of a 2B-FLX product "
                       " granule (using surface type from AUX-MET + AUX-SAT)")
    arg_parser = argparse.ArgumentParser(description=arg_description)
    arg_parser.add_argument("-s", "--use_saved_exe", action="store_true",
                            help="Use the MatLab Runtime to execute routines "
                                 "stored in a Matlab executable file.")
    arg_parser.add_argument("-i", "--interactive-no_display",
                            dest="interactive", action="store_true",
                            help="Set environment, then run MatLab in "
                               "interactive mode (but with no fancy display).") 

    args = arg_parser.parse_args()

    # Check for argument sanity:
    if args.use_saved_exe and args.interactive:
        print ("ERROR: Arguments -s and -i cannot be used together.")
        sys.exit(1)

    # Run driver:
    main(args.use_saved_exe, args.interactive)
