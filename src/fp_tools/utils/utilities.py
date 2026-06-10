#!/usr/bin/env python

"""Utilities for argparse, logging, file handling, and multiprocessing."""

import os
import logging
import sys
import re
import argparse
import collections
import time
import textwrap
import string
import numpy as np
import multiprocessing as mp
import copy
from difflib import SequenceMatcher
import traceback

import pyBigWig
from fp_tools.utils.logger import *


PROGRESS_LOG_STEP = 10


def progress_log_percent(done, total, previous_percent=None, step=PROGRESS_LOG_STEP):
    """Return the next coarse progress percent to log, or None if unchanged."""

    if total <= 0:
        return None

    percent = int(round(100.0 * done / float(total)))
    percent = max(0, min(100, percent))
    if percent >= 100:
        milestone = 100
    else:
        milestone = int(percent / step) * step

    if previous_percent is None or milestone > previous_percent:
        return milestone
    return None


def show_worker_progress(verbosity, total_items, is_tty=None):
    """Only show per-worker progress bars for explicit high-verbosity TTY runs."""

    if is_tty is None:
        is_tty = sys.stderr.isatty()
    return bool(verbosity >= 5 and total_items > 1 and is_tty)


#-------------------------------------------------------------------------------------------#
#----------------------------------- Multiprocessing ---------------------------------------#
#-------------------------------------------------------------------------------------------#
def check_cores(given_cores, logger):
    """Resolve the effective core count for local execution.

    Positive user-supplied values are respected and capped at the local CPU
    count. Invalid values fall back to all available cores.
    """

    available_cores = mp.cpu_count()

    if given_cores is None or given_cores < 1:
        logger.warning(
            "Invalid '--cores' value {0}; using all available cores ({1}).".format(
                given_cores, available_cores
            )
        )
        return available_cores

    if given_cores > available_cores:
        logger.warning(
            "Requested {0} cores, but only {1} are available; using {1}.".format(
                given_cores, available_cores
            )
        )
        return available_cores

    return given_cores


def run_parallel(FUNC, input_chunks, arguments, n_cores, logger, progress_label="Progress:"):
    """
    Run FUNC over input_chunks in parallel, preserving submit order in results.
    - FUNC is called as: FUNC(chunk, *arguments)
    - Progress prints like the original ("Progress: XX%"), or uses tqdm when a TTY is available.
    - Returns: list of FUNC outputs in the same order as input_chunks
    """
    no_chunks = len(input_chunks)
    if no_chunks == 0:
        return []

    # Try tqdm (pretty) only when we're in a TTY; otherwise fall back to the original log lines
    try:
        from tqdm import tqdm
        use_tqdm = (sys.stderr.isatty() or sys.stdout.isatty()) and no_chunks > 0
    except Exception:
        tqdm = None
        use_tqdm = False

    # -------- Parallel path --------
    if (n_cores or 1) > 1:
        pool = mp.Pool(processes=n_cores)
        try:
            task_list = [pool.apply_async(FUNC, args=[chunk] + list(arguments)) for chunk in input_chunks]
            pool.close()  # no more tasks submitted

            if use_tqdm:
                # tqdm to stderr to avoid mixing with logger stdout
                bar = tqdm(total=no_chunks, desc=progress_label.rstrip(' :'), unit="task", leave=True, file=sys.stderr)
                done_prev = -1
                try:
                    while True:
                        done = sum(1 for t in task_list if t.ready())
                        if done != done_prev:
                            bar.n = done
                            bar.refresh()
                            done_prev = done
                        if done >= no_chunks:
                            break
                        time.sleep(0.5)  # keep the original pacing
                finally:
                    bar.close()
                    # ensure the next stdout log starts on a fresh line
                    try:
                        logger.info("")
                    except Exception:
                        pass
            else:
                # Coarse log milestones avoid one line per percent in batch logs.
                last_done = -1
                last_pct = None
                while True:
                    done = sum(1 for t in task_list if t.ready())
                    if done != last_done:
                        pct = progress_log_percent(done, no_chunks, last_pct)
                        if pct is not None:
                            logger.info(f"{progress_label} {pct}%")
                            last_pct = pct
                        last_done = done
                    if done >= no_chunks:
                        break
                    time.sleep(0.5)

            # After completion, join the pool and collect results in submit order
            pool.join()
            output_list = [task.get() for task in task_list]
            return output_list

        finally:
            # Be safe if an exception occurred above
            try:
                pool.terminate()
            except Exception:
                pass
            try:
                pool.join()
            except Exception:
                pass

    # -------- Single-core path (preserve original semantics) --------
    output_list = []
    if use_tqdm:
        from tqdm import tqdm  # already imported above; local for clarity
        for idx, chunk in enumerate(tqdm(input_chunks, desc=progress_label.rstrip(' :'), unit="task", leave=True, file=sys.stderr)):
            # keep logger behavior compatible: print an empty line after tqdm finishes
            if idx == no_chunks - 1:
                pass
            output_list.append(FUNC(chunk, *arguments))
        try:
            logger.info("")
        except Exception:
            pass
    else:
        last_pct = None
        for count, input_chunk in enumerate(input_chunks):
            pct = progress_log_percent(count, no_chunks, last_pct)
            if pct is not None:
                logger.info("Progress: {0}%".format(pct))
                last_pct = pct
            output_list.append(FUNC(input_chunk, *arguments))
        # ensure we also log the final 100% like before
        pct = progress_log_percent(no_chunks, no_chunks, last_pct)
        if pct is not None:
            logger.info("Progress: 100%")

    return output_list


def file_writer(q, key_file_dict, args):
    """ File-writer per key -> to value file """

    #Open handles for all files (but only once per file!)
    file2handle = {}
    for fil in set(key_file_dict.values()):
        try:
            file2handle[fil] = open(fil, "w")
        except Exception as e:
            print("Error opening file {0} in file_writer. Exception was: '{1}'".format(fil, e))
            raise e
            return(1)

    #Assign handles to keys
    handles = {}
    for key in key_file_dict:
        handles[key] = file2handle[key_file_dict[key]]

    #Fetching string content from queue
    while True:
        try:
            (key, content) = q.get()
            if key == None:
                break

            handles[key].write(content)

        except Exception as e:
            import sys, traceback
            print('Problem in file_writer: ')
            print(e)
            raise e
            break

    #Got all regions in queue, close files
    for fil, handle in file2handle.items():
        handle.close()

    return(0)	#success


def bigwig_writer(q, key_file_dict, header, regions, args):
    """ Handle queue to write bigwig, args contain extra info such as verbosity and log_q """

    #todo: check if args.log_q exists
    logger = FpToolsLogger("", args.verbosity, args.log_q)	#separate core, needs mp logger through queue
    logger.debug("Opened bigwig writer process for {0}".format(key_file_dict))
    logger.debug("Header: {0}".format(header))

    handles = {}
    for key in key_file_dict:
        logger.debug("Opening file {0} for writing".format(key_file_dict[key]))
        try:
            handles[key] = pyBigWig.open(key_file_dict[key], "w")
            handles[key].addHeader(header)

        except Exception as e:
            logger.error("Error opening file {0} in bigwig_writer. Exception was: '{1}'".format(key_file_dict[key], e))
            raise e

    #Correct order of chromosomes as given in header
    contig_list = [tup[0] for tup in header]
    order_dict = dict(zip(contig_list, range(len(contig_list))))

    #Establish order of regions to be writteninput regions
    region_tups = [(region.chrom, region.start, region.end) for region in regions]
    sorted_region_tups = sorted(region_tups, key=lambda tup: (order_dict[tup[0]], tup[1]))			#sort to same order as bigwig header
    n_regions = len(region_tups)

    #Fetching content from queue
    logger.debug("Fetching content from queue")

    i_to_write = {key:0 for key in handles}			#index of next region to write
    ready_to_write = {key:{} for key in handles}	#key:dict; dict is region-tup:signal array
    while True:

        try:
            (key, region, signal) = q.get()	#key is the bigwig key (e.g. bias:forward), region is a tup of (chr, start, end)
            logger.spam("Received signal {0} for region {1}".format(key, region))

            if key == None:	#none is only inserted once all regions have been sent
                for akey in i_to_write:
                    if i_to_write[akey] != n_regions:	#i_to_write is written index + 1
                        logger.error("Wrote {0} regions but there are {1} in total".format(i_to_write[akey], len(region_tups)))
                        logger.error("Ready_to_write[{0}]: {1}".format(akey, len(ready_to_write[akey])))
                        sys.exit()
                break

            #Save key:region:signal to ready_to_write
            ready_to_write[key][region] = signal

            writing_progress = Progress(n_regions, logger, prefix="Writing progress", round=0)

            #Check if next-to-write region was done
            for key in handles:
                logger.spam("Key: {0}. Index to write: {1}".format(key, i_to_write[key]))

                #Only deal with writing if there are still regions to write for this handle
                if i_to_write[key] < n_regions: 	#if i_to_write == 2 and n_regions == 2, the two regions (idx 0+1) have already been written
                    next_region = sorted_region_tups[i_to_write[key]]	#this is the region to be written next for this key

                    #If results are in; write wanted entry to bigwig
                    while next_region in ready_to_write[key]: 	#When true: Keep writing when the next region is available
                        chrom = next_region[0]
                        signal = ready_to_write[key][next_region]
                        included = signal.nonzero()[0]
                        positions = np.arange(next_region[1],next_region[2])		#start-end	(including end)
                        pos = positions[included].tolist()
                        val = signal[included].tolist()

                        if len(pos) > 0:
                            try:
                                handles[key].addEntries(chrom, pos, values=val, span=1)
                            except Exception as e:
                                logger.error("Error writing key: {0}, region: {1} to bigwig. Exception was: '{2}'".format(key, next_region, e))
                                logger.debug("Chrom: {0}".format(chrom))
                                logger.debug("Positions: {0}".format(pos))
                                logger.debug("Values: {0}".format(val))
                                raise e
                        logger.spam("Wrote signal {0} from region {1}".format(key, next_region))

                        #Free up memory in dict
                        ready_to_write[key][next_region] = None

                        #Check whether this was the last region
                        if i_to_write[key] == n_regions - 1: #i_to_write is the last idx in regions; all sites were written
                            logger.info("Closing {0} (this might take some time)".format(key_file_dict[key]))
                            handles[key].close()
                            i_to_write[key] += 1
                            next_region = None 	#exit the while loop
                        else:
                            i_to_write[key] += 1
                            next_region = sorted_region_tups[i_to_write[key]]	#this is the region to be written next for this key

                        #Writing progress
                        #progress = sum([i_to_write[key] for key in handles])
                        #writing_progress.write(progress)

        except Exception as e:
            logger.error("Problem in bigwig_writer. Exception was: '{0}'".format(e))
            traceback.print_tb(e.__traceback__)
            raise e
            return(1) #Return with error

    return(0)	#everything went well


def monitor_progress(task_list, logger, prefix="Progress"):
    """
    Show progress for a list of multiprocessing AsyncResult tasks.
    Uses tqdm (to stderr) when available and attached to a TTY, else logs lines.
    """
    total = len(task_list)

    try:
        from tqdm import tqdm
        use_tty = sys.stderr.isatty()
        use_tqdm = bool(tqdm) and use_tty and total > 0
    except Exception:
        tqdm = None
        use_tqdm = False

    if use_tqdm:
        bar = tqdm(total=total, desc=prefix.rstrip(' :'), unit="task", leave=True, file=sys.stderr)
        done_prev = -1
        try:
            while True:
                done = sum(1 for t in task_list if t.ready())
                if done != done_prev:
                    bar.n = done
                    bar.refresh()
                    done_prev = done
                if done >= total:
                    break
                time.sleep(0.2)
        finally:
            bar.close()
            # ensure the next stdout log starts on a fresh line
            try:
                logger.info("")  # harmless empty line
            except Exception:
                pass
    else:
        last_pct = None
        while True:
            done = sum(1 for t in task_list if t.ready())
            pct = progress_log_percent(done, total, last_pct)
            if pct is not None:
                logger.info(f"{prefix} {pct}%")
                last_pct = pct
            if done >= total:
                break
            time.sleep(0.2)
#-------------------------------------------------------------------------------------------#
#------------------------------------- Argparser -------------------------------------------#
#-------------------------------------------------------------------------------------------#

def restricted_float(f, f_min, f_max):
    f = float(f)
    if f < f_min or f > f_max:
        raise argparse.ArgumentTypeError("{0} not in range [0.0, 1.0]".format(f))
    return f

def restricted_int(integer, i_min, i_max):
    integer = float(integer)
    if integer < i_min or integer > i_max:
        raise

def format_help_description(name, description, width=90):
    """ Format description of command line tool --help description """

    formatted = "" #initialize

    #Calculate needed whitespace in comparison to header
    header = "fp-tools {0}".format(name)
    ws = int((width - len(header))/2.0)

    formatted += "_"*width + "\n"*2
    formatted += "{0}{1}{0}\n".format(" "*ws, header)
    formatted += "_"*width + "\n"*2

    for segment in description.split("\n"):
        formatted += "\n".join(textwrap.wrap(segment, width)) + "\n"	#Split description on space

    if description != "":
        formatted += "\n" + "-"*width + "\n"

    return(formatted)


def check_required(args, required):
    """ Checks required keys in input args """

    for arg in required:
        if getattr(args, arg) == None:
            sys.exit("ERROR: Missing argument --{0}".format(arg))

def add_underscore_options(parser):

    for group in parser._action_groups:
        group_actions = group._group_actions

        if len(group_actions) > 0:
            for option in group_actions:
                opt_string = option.option_strings[-1]
                opt_string_fmt = re.sub(r'^\-*', "", opt_string)

                #Add backwards compatibility of options with -/_
                if "-" in opt_string_fmt:
                    new_opt_string = "--" + opt_string_fmt.replace("-", "_")

                    #Get keys for the new option
                    keep = ["nargs", "const", "default", "type", "choices", "required", "metavar"]
                    new_option_dict = {key: option.__dict__[key] for key in keep}
                    new_option_dict["nargs"] = "?" if new_option_dict["nargs"] == 0 else new_option_dict["nargs"]

                    parser.add_argument(new_opt_string, help=argparse.SUPPRESS, **new_option_dict)

    return(parser)

#-------------------------------------------------------------------------------------------#
#---------------------------------------- Misc ---------------------------------------------#
#-------------------------------------------------------------------------------------------#

def num(s):
    try:
        return int(s)
    except ValueError:
        return float(s)

class Progress:
    """ Class for writing out progress of processes such as multiprocessing """

    def __init__(self, total_elements, logger, prefix="Progress", round=0):

        self.total_elements = total_elements
        self.logger = logger
        self.prefix = prefix
        self.round = round
        self.last_written = None

    def write(self, progress):
        """ Write out progress if it was not already written """

        if self.round == 0:
            percent_to_write = progress_log_percent(progress, self.total_elements, self.last_written)
        else:
            raw_percent = progress / self.total_elements * 100
            percent_to_write = round(raw_percent, self.round)
            if percent_to_write == self.last_written:
                percent_to_write = None

        #Only write if this level has not already been written
        if percent_to_write is not None:
            self.logger.info("{0}: {1}%".format(self.prefix, percent_to_write))
            self.last_written = percent_to_write


def flatten_list(lst):

    for element in lst:
        if hasattr(element, "__iter__") and not isinstance(element, (str, bytes)):
            yield from flatten_list(element)
        else:
            yield element

def expand_dirs(list_of_paths):
    """ Expands a list of files and dirs to a list of all files within dirs """

    all_files = []
    for path in list_of_paths:
        if os.path.isdir(path):
            files = os.listdir(path)
            all_files.extend([os.path.join(path, f) for f in files])
        else:
            all_files.append(path)

    return(all_files)

def check_files(lst_of_files, action="r", logger=FpToolsLogger()):

    flat_lst = flatten_list(lst_of_files)
    for fil in flat_lst:
        if fil != None:
            logger.debug("Checking " + fil)
            if action == "r":
                if os.path.exists(fil):
                    try:
                        with open(fil) as f:
                            pass
                    except:
                        sys.exit("ERROR: Could not open file \"{0}\" for reading".format(fil))
                else:
                    sys.exit("ERROR: File \"{0}\" does not exists".format(fil))

            elif action == "w":
                if os.path.exists(fil):
                    try:
                        with open(fil, "w") as f:
                            pass
                    except:
                        sys.exit("ERROR: Could not open file \"{0}\" for writing. Please check that you do not have the file open and that you have write permission.".format(fil))

def make_directory(directory):
    if not os.path.isfile(directory) and not os.path.exists(directory):
        os.makedirs(directory)


def merge_dicts(dicts):
    """ Merge recursive keys and values for list of dicts into one dict. Values are added numerically / lists are extended / numpy arrays are added"""

    def merger(dct, dct_to_add):
        """ Merge recursive keys and values of dct_to_add into dct. Values are added numerically / lists are extended / numpy arrays are added
        No return - dct is changed in place """

        for k, v in dct_to_add.items():

                #If k is a dict, go one level down
                if (k in dct and isinstance(dct[k], dict)):
                    merger(dct[k], dct_to_add[k])
                else:
                    if not k in dct:
                        dct[k] = dct_to_add[k]
                    else:
                        dct[k] += dct_to_add[k]

    #Initialize with the first dict in list
    out_dict = copy.deepcopy(dicts[0])
    for dct in dicts[1:]:
        merger(out_dict, dct)

    return(out_dict)

def filafy(astring):
    """ Make string into accepted filename """

    valid_chars = "-_.%s%s" % (string.ascii_letters, string.digits)
    filename = ''.join(char for char in astring if char in valid_chars)
    return(filename)


def get_closest(value, arr): 
    """Find element the element in arr which is closest to value """

    idx = (np.abs(arr-value)).argmin()
    return(arr[idx])


#-------------------------------------------------------------------------------------------#
#-------------------------------------- Matching -------------------------------------------#
#-------------------------------------------------------------------------------------------#

def common_prefix(strings):
    """ Find the longest string that is a prefix of all the strings. Used in PlotChanges to find TF names from list """
    if not strings:
        return ''
    prefix = strings[0]
    for s in strings:
        if len(s) < len(prefix):
            prefix = prefix[:len(s)]
        if not prefix:
            return ''
        for i in range(len(prefix)):
            if prefix[i] != s[i]:
                prefix = prefix[:i]
                break
    return prefix


def match_lists(lofl): # list of lists
    """ Find matches between list1 and list2 (output will be the length of list1 with one or more matches per element)."""

    #Remove common prefixes/suffixes per list
    prefixes = []
    suffixes = []
    for i, lst in enumerate(lofl):
        if len(lst) > 1:	#make sure to only compare between more than one element (otherwise it sets the whole element as prefix)
            prefix = common_prefix(lst)
            suffix = common_prefix([element[::-1] for element in lst])[::-1]
        else:
            prefix = ""
            suffix = ""

        lofl[i] = [re.sub("^" + prefix, "", re.sub(suffix + "$", "", element)) for element in lst]
        prefixes.append(prefix)
        suffixes.append(suffix)

    #Initialize
    matches = [[[] for _ in lofl[0]] for _ in range(len(lofl)-1)]	#list of lists (len(lofl) * len(lofl[0]) * list for collecting results)

    #Find best matches to each element in col:
    for col in range(len(lofl)-1):

        lst1 = lofl[0]  	#[row for row in lofl[0]]	#Always match to first list
        lst2 = lofl[col+1] 	# [row for row in lofl[col+1]] 	#compare to column

        for i, element1 in enumerate(lst1):

            local_match_scores = []
            global_match_scores = []
            for element2 in lst2:
                match_tups = SequenceMatcher(None, element1.lower(), element2.lower()).get_matching_blocks()
                match_sum = sum([tup.size for tup in match_tups if tup.size > 1])

                local_match_scores.append(match_sum / float(len(element1)))		#local score shows how well element1 fits locally in element2
                global_match_scores.append(match_sum / float(len(element1) + len(element2)))	#takes into account length of both element1 and 2

            #Best local match scores
            best_score = max(local_match_scores)
            best_idx = [idx for idx, score in enumerate(local_match_scores) if score == best_score]

            #Sort best_idx by global match scores
            best_global_match_scores = [global_match_scores[idx] for idx in best_idx]
            best_idx = [idx for _,idx in sorted(zip(best_global_match_scores, best_idx), key=lambda pair: pair[0], reverse=True)]	#largest global match scores first
            best_elements = [lst2[idx] for idx in best_idx]		#elements in lst2 that fit to element1
            matches[col][i].extend(best_elements)

    #Map back to full names
    for col in range(1, len(lofl)): 	# number of lists to get matches from
        for row in range(len(lofl[0])): 	#elements in lofl[0]
            matches[col-1][row] = [prefixes[col] + match + suffixes[col] for match in matches[col-1][row]]

    return(matches)
