from fp_tools.parsers import add_aggregate_arguments
from fp_tools.tools.plot_aggregate import run_aggregate


def main():
    import argparse
    parser = add_aggregate_arguments(argparse.ArgumentParser())
    args = parser.parse_args()
    run_aggregate(args)
