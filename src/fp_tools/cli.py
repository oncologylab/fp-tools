from fp_tools.parsers import add_atacorrect_arguments
from fp_tools.tools.atacorrect import run_atacorrect

def main():
    import argparse
    parser = add_atacorrect_arguments(argparse.ArgumentParser())
    args = parser.parse_args()
    run_atacorrect(args)
