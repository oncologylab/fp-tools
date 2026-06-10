from fp_tools.parsers import add_scorebigwig_arguments
from fp_tools.tools.score_bigwig import run_scorebigwig

def main():
    import argparse
    parser = add_scorebigwig_arguments(argparse.ArgumentParser())
    args = parser.parse_args()
    run_scorebigwig(args)
