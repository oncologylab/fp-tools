# src/fp_tools/utils/logger.py
#!/usr/bin/env python

"""Shared logging helpers for fp-tools command-line workflows."""

import sys
import os
import re
from datetime import datetime
import logging
import logging.handlers
import multiprocessing as mp
import time
from importlib.metadata import PackageNotFoundError, version as package_version

DISPLAY_NAME = "fp-tools"
DIST_NAME = "fp-tools-bio"

def add_logger_args(args):
    """Add verbosity option to argparse parser"""
    args.add_argument(
        '--verbosity', metavar="<int>",
        help="Level of output logging (0: silent, 1: errors/warnings, 2: info, 3: stats, 4: debug, 5: spam) (default: 3)",
        choices=[0,1,2,3,4,5], default=3, type=int
    )
    return args

class FpToolsLogger(logging.Logger):
    """Logger with convenience helpers used across fp-tools commands."""

    logger_levels = {
        0: 0,
        1: logging.WARNING,                         # warnings + errors
        2: logging.INFO,                            # info
        3: int((logging.INFO + logging.DEBUG) / 2), # stats (between info and debug)
        4: logging.DEBUG,                           # debug
        5: logging.DEBUG - 5                        # spam-level debug
    }

    def __init__(self, tool_name=DISPLAY_NAME, level=3, queue=None):
        self.tool_name = tool_name
        super().__init__(self.tool_name)

        if level == 0:
            self.disabled = True

        # custom levels
        comment_level = FpToolsLogger.logger_levels[1] + 1
        logging.addLevelName(comment_level, "comment")
        setattr(self, "comment", lambda *args: self.log(comment_level, *args))

        stats_level = FpToolsLogger.logger_levels[3]
        logging.addLevelName(stats_level, "STATS")
        setattr(self, "stats", lambda *args: self.log(stats_level, *args))

        spam_level = FpToolsLogger.logger_levels[5]
        logging.addLevelName(spam_level, "SPAM")
        setattr(self, "spam", lambda *args: self.log(spam_level, *args))

        self.level = FpToolsLogger.logger_levels[level]
        self.formatter = FpToolsFormatter()
        self.setLevel(self.level)

        if queue is None:
            con = logging.StreamHandler(sys.stdout)
            con.setLevel(self.level)
            con.setFormatter(self.formatter)
            self.addHandler(con)
        else:
            h = logging.handlers.QueueHandler(queue)
            self.handlers = []
            self.addHandler(h)

        self.begin_time = datetime.now()
        self.end_time = None
        self.total_time = None
        self.queue = None
        self.listener = None

    def _version_string(self):
        ver = None
        try:
            ver = package_version(DIST_NAME)
        except PackageNotFoundError:
            pass
        try:
            from fp_tools import __version__ as FP_VERSION
            ver = FP_VERSION or ver
        except Exception:
            pass
        return ver or "unknown"

    def begin(self):
        """Write header lines for the run"""
        version = self._version_string()
        # executable name, fall back to the logical tool name
        import os, sys as _sys
        prog = os.path.basename(_sys.argv[0]) or self.tool_name or "fp-tools"
        # header
        self.comment(f"# {DISPLAY_NAME} {version} {self.tool_name} (run started {self.begin_time})")
        self.comment(f"# Working directory: {os.getcwd()}")
        self.comment(f"# Command line call: {prog} " + " ".join(_sys.argv[1:]) + "\n")

    def stop(self):
        self.end_time = datetime.now()
        self.total_time = self.end_time - self.begin_time

    def end(self):
        self.end_time = datetime.now()
        self.total_time = self.end_time - self.begin_time
        self.comment("")
        self.info(f"Finished {self.tool_name} run (total time elapsed: {self.total_time})")

    def start_logger_queue(self):
        self.debug("Starting logger queue for multiprocessing")
        self.queue = mp.Manager().Queue()
        self.listener = mp.Process(target=self.main_logger_process)
        self.listener.start()

    def stop_logger_queue(self):
        self.debug("Waiting for listener to finish")
        if self.queue is not None:
            self.queue.put(None)
        if self.listener is not None:
            while self.listener.exitcode != 0:
                self.debug(f"Listener exitcode is: {self.listener.exitcode}. Waiting for exitcode = 0.")
                time.sleep(0.1)
            self.debug("Joining listener")
            self.listener.join()

    def main_logger_process(self):
        self.debug("Started main logger process")
        while True:
            try:
                record = self.queue.get()
                if record is None:
                    break
                self.handle(record)
            except EOFError:
                self.error("Multiprocessing logger lost connection to queue - probably due to an error raised from a child process.")
                break
        return 1

    # ---- overview helpers used by command drivers ----
    def arguments_overview(self, parser, args):
        """Print a simple overview of parsed CLI arguments."""
        content = ""
        content += "# ----- Input parameters -----\n"
        for group in parser._action_groups:
            group_actions = group._group_actions
            if len(group_actions) > 0:
                for option in group_actions:
                    if option.help != "==SUPPRESS==":
                        name = option.dest
                        attr = getattr(args, name, None)
                        content += f"# {name}:\t{attr}\n"
        self.comment(content + "\n")

    def output_files(self, outfiles):
        """Print the list of output files."""
        self.comment("# ----- Output files -----")
        for outf in outfiles:
            if outf is not None:
                self.comment(f"# {outf}")
        self.comment("\n")


class FpToolsFormatter(logging.Formatter):
    """Formatter used by FpToolsLogger"""
    default_fmt = logging.Formatter(
        "%(asctime)s (%(process)d) [%(levelname)s]\t%(message)s", "%Y-%m-%d %H:%M:%S"
    )
    comment_fmt = logging.Formatter("%(message)s")

    def format(self, record):
        if record.levelname == "comment":
            return self.comment_fmt.format(record)
        elif record.levelno != 0:
            return self.default_fmt.format(record)
        else:
            return ""
