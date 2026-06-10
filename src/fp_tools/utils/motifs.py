#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Utilities for reading, writing, plotting and scanning DNA motifs.

The module provides package-local motif containers, plotting helpers,
format readers, and MOODS-based scanning utilities.
"""

from __future__ import annotations

import base64
import io
import math
import os
import re
import sys
import copy
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

# --- Optional deps ---------------------------------------------------------
try:  # plotting style; optional
    import seaborn as sns  # type: ignore
except Exception:  # pragma: no cover
    sns = None

try:  # logos
    import logomaker  # type: ignore
except Exception as e:  # pragma: no cover
    raise ImportError(
        "logomaker is required for motif logos. Install with: `poetry add logomaker`"
    ) from e

try:  # BioPython motifs reader/writer
    from Bio import motifs as bio_motifs  # type: ignore
except Exception as e:
    raise ImportError(
        "biopython is required for reading/writing motifs. Install with: `poetry add biopython`"
    ) from e

try:  # fast scanning
    import MOODS.scan  # type: ignore
    import MOODS.tools  # type: ignore
    import MOODS.parsers  # type: ignore
except Exception as e:  # pragma: no cover
    raise ImportError(
        "MOODS is required for motif scanning. Install with: `poetry add MOODS-python==1.9.4.1`"
    ) from e

# --- Internal (fp_tools namespace) ----------------------------------------
from fp_tools.utils.regions import OneRegion, RegionList

# We try to import helpers from your utilities, but also provide local
# fallbacks so this module works out-of-the-box.
try:
    from fp_tools.utils.utilities import filafy as _filafy, num as _num  # type: ignore
except Exception:  # pragma: no cover
    _filafy = None
    _num = None


def filafy(s: str) -> str:
    """Make a filesystem-friendly string (fallback implementation)."""
    if _filafy is not None:
        return _filafy(s)
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", str(s)).strip("_")
    return s or "None"


def num(x):
    """Best-effort numeric cast (fallback implementation)."""
    if _num is not None:
        return _num(x)
    try:
        i = int(x)
        return i
    except Exception:
        try:
            return float(x)
        except Exception:
            return x


# --------------------------------------------------------------------------
# helpers

def float_to_int(afloat):
    """Converts integer-like floats (e.g. 1.0000) to ints; otherwise returns original."""
    parts = str(afloat).split(".")
    if len(parts) == 1:
        return afloat
    if len(parts) == 2 and float(parts[1]) == 0:
        return int(float(afloat))
    return afloat


def is_symmetric(matrix: np.ndarray) -> bool:
    """Check if a matrix is symmetric around the diagonal."""
    s = matrix.shape
    if s[0] != s[1]:
        return False
    return np.allclose(matrix, matrix.T, equal_nan=True)


# --------------------------------------------------------------------------
# Core data structures

class MotifList(list):
    """A list of :class:`OneMotif` objects with utilities for scanning/clustering."""

    # these are filled by setup_moods_scanner
    forward: Dict[str, list] = {"names": [], "matrices": [], "thresholds": []}
    reverse: Dict[str, list] = {"names": [], "matrices": [], "thresholds": []}

    moods_scanner_forward = None
    moods_scanner_reverse = None

    def __init__(self, lst: Optional[List["OneMotif"]] = None):
        super().__init__(iter(lst or []))
        self.set_background()

    def __str__(self) -> str:  # pragma: no cover
        return "\n".join([str(onemotif) for onemotif in self])

    # ---------------- I/O --------------------------------------------------
    def from_file(self, path: str) -> "MotifList":
        """Read motifs from a file into this MotifList.

        Supports "meme", "jaspar", "pfm" and "transfac" formats.
        """
        content = open(path).read()
        file_format = get_motif_format(content)

        if file_format == "meme":
            bg_flag = False
            proba_flag = False
            bases = None
            strands = None
            bg = None
            lines = content.split("\n") + [""]
            for line in lines:
                if line.startswith("ALPHABET="):
                    bases = list(line.replace("ALPHABET=", "").strip())
                elif line.startswith("strands"):
                    strands = line.replace("strands: ", "")
                elif line.startswith("Background letter frequencies"):
                    bg_flag = True
                elif bg_flag:
                    bg_flag = False
                    bg_and_freq = re.split(r"(?<=\d)\s+", line.strip())
                    bg_dict = {k: float(v) for k, v in [el.split(" ") for el in bg_and_freq]}
                    if bases is None:
                        bases = sorted(bg_dict.keys())
                    bg = np.array([bg_dict[base] for base in bases])
                elif line.startswith("MOTIF"):
                    if proba_flag:
                        self[-1].pfm = np.array(probability_matrix).T
                        count_matrix = (self[-1].pfm * self[-1].n)
                        count_matrix = np.round(count_matrix).astype(int).tolist()
                        self[-1].set_counts(count_matrix)
                    probability_matrix = []  # type: ignore
                    self.append(OneMotif(motifid=""))
                    cols = line.split()
                    if len(cols) > 2:
                        motif_id, name = cols[1], cols[2]
                    else:
                        motif_id, name = cols[1], ""
                    self[-1].id = motif_id
                    self[-1].name = name
                    if bases is not None:
                        self[-1].bases = bases
                    if strands is not None:
                        self[-1].strands = strands
                    if bg is not None:
                        self[-1].bg = bg
                elif line.startswith("letter-probability matrix"):
                    proba_flag = True
                    key_value_string = re.sub("letter-probability matrix:\\s*", "", line.rstrip())
                    if len(key_value_string) > 0:
                        key_value_split = re.split(r"(?<!=)\s+", key_value_string)
                        key_value_lists = [re.split(r"=\s*", pair) for pair in key_value_split]
                        key_value_dict = {pair[0]: pair[1] for pair in key_value_lists}
                        info = key_value_dict
                        self[-1].info = info
                        if "nsites" in info:
                            self[-1].n = int(float(info["nsites"]))
                elif proba_flag:
                    if re.match(r"^\s*(?![-\s]+)([\d\-\.\se\-]+?)$", line):
                        columns = list(map(float, line.split()))
                        if not len(columns) == len(self[-1].bases):
                            sys.exit(
                                "Error when reading probability matrix from {0}! Expected {1} columns found {2}!".format(
                                    path, len(self[-1].bases), len(columns)
                                )
                            )
                        probability_matrix.append(columns)
                    else:
                        proba_flag = False
                        self[-1].pfm = np.array(probability_matrix).T
                        count_matrix = (self[-1].pfm * self[-1].n)
                        count_matrix = np.round(count_matrix).astype(int).tolist()
                        self[-1].set_counts(count_matrix)

        elif file_format in ["jaspar", "pfm"]:
            with open(path) as f:
                for line in f:
                    if line.startswith(">"):
                        if len(self) > 0:
                            self[-1].set_counts(count_matrix)  # type: ignore
                        cols = line.split()
                        if len(cols) > 1:
                            motif_id, name = cols[0].replace(">", ""), cols[1]
                        else:
                            motif_id, name = cols[0].replace(">", ""), ""
                        self.append(OneMotif(name=name, motifid=motif_id))
                        count_matrix = []  # type: ignore
                    else:
                        if len(self) == 0:
                            raise ValueError(
                                f"Error when reading motifs from {path}! No motif header found before first motif."
                            )
                        if len(line.rstrip()) == 0:
                            continue
                        values = (
                            re.sub(r"^[ACGT\s\[]+", "", line.strip()).replace("]", "").split()
                        )
                        counts = list(map(float, values))
                        count_matrix.append(counts)  # type: ignore
                if len(self) > 0:
                    self[-1].set_counts(count_matrix)  # type: ignore
                for motif in self:
                    motif.get_pfm()

        elif file_format == "transfac":
            with open(path) as f:
                for m in bio_motifs.parse(f, file_format, strict=False):
                    self.append(
                        OneMotif(
                            motifid=m.get("AC", ""),
                            name=m.get("ID", ""),
                            counts=[m.counts[base] for base in ["A", "C", "G", "T"]],
                        )
                    )
                    self[-1].biomotifs_obj = m  # type: ignore[attr-defined]
        else:
            sys.exit(f"Error when reading motifs from {path}! Unsupported file format: {file_format}")

        # validate widths
        for motif in self:
            if "w" in motif.info:
                w = int(motif.info["w"])  # type: ignore[arg-type]
                l = len(motif.counts[0])
                if w != l:
                    sys.exit(
                        "Error reading motif '{3}' from {0}! 'w' given in header ({1}) does not match motif length ({2})".format(
                            path, w, l, motif.id
                        )
                    )

        for motif in self:
            if motif.length is None:
                sys.exit(
                    f"ERROR: No matrix could be read for motif '{motif.id} {motif.name}' - please check the input file."
                )

        # intify
        for motif in self:
            for r in range(4):
                for c in range(motif.length):
                    motif.counts[r][c] = float_to_int(motif.counts[r][c])
        return self

    def to_file(self, path: str, fmt: str = "pfm") -> "MotifList":
        with open(path, "w") as f:
            f.write(self.as_string(fmt))
        return self

    def as_string(self, output_format: str = "pfm") -> str:
        out = []
        header = True
        for motif in self:
            out.append(motif.as_string(output_format=output_format, header=header))
            header = False
        return "".join(out)

    # --------------- background ------------------------------------------
    def get_background(self) -> Optional[np.ndarray]:
        if len(self) == 0:
            return None
        global_bg = np.array([0.0] * len(self[0].bg))
        total_n = 0
        for motif in self:
            global_bg += motif.bg * motif.n
            total_n += motif.n
        return global_bg / total_n

    def set_background(self) -> None:
        self.bg = self.get_background() if len(self) > 0 else np.array([0.25] * 4)

    # --------------- scanning (MOODS) ------------------------------------
    def setup_moods_scanner(self, strand: str = ".") -> None:
        for motif in self:
            if len(motif.prefix) <= 0 or motif.threshold is None:
                raise Exception(
                    "Missing prefix and/or threshold! Run motif.set_prefix() and motif.get_threshold() first."
                )
        if strand in ["+", "."]:
            self.moods_scanner_forward, self.forward_parameters = self.__init_scanner(strand="+")
        if strand in ["-", "."]:
            self.moods_scanner_reverse, self.reverse_parameters = self.__init_scanner(strand="-")

    def __init_scanner(self, strand: str = "+"):
        if strand == "+":
            motifs = self
        else:
            motifs = self.get_reverse()
            motifs.set_background()
        for motif in motifs:
            if motif.pssm is None:
                motif.get_pssm()
        parameters = {"names": [], "matrices": [], "thresholds": []}
        for motif in motifs:
            parameters["names"].append(motif.prefix)
            parameters["matrices"].append(motif.pssm)
            parameters["thresholds"].append(motif.threshold)
        scanner = MOODS.scan.Scanner(7)
        scanner.set_motifs(parameters["matrices"], motifs.bg, parameters["thresholds"])  # type: ignore[attr-defined]
        return (scanner, parameters)

    def scan_sequence(self, seq: str, region: OneRegion, strand: str = ".") -> RegionList:
        sites = RegionList()
        if strand in ["+", "."] and self.moods_scanner_forward is None:
            self.setup_moods_scanner("+")
        if strand in ["-", "."] and self.moods_scanner_reverse is None:
            self.setup_moods_scanner("-")
        if strand in ["+", "."]:
            sites += self.__stranded_scan(seq=seq, region=region, strand="+")
        if strand in ["-", "."]:
            sites += self.__stranded_scan(seq=seq, region=region, strand="-")
        return sites

    def __stranded_scan(self, seq: str, region: OneRegion, strand: str = "+") -> RegionList:
        sites = RegionList()
        if strand == "+":
            scanner = self.moods_scanner_forward
            parameters = self.forward_parameters
        else:
            scanner = self.moods_scanner_reverse
            parameters = self.reverse_parameters
        results = scanner.scan(seq)  # type: ignore[union-attr]
        for (matrix, name, result) in zip(parameters["matrices"], parameters["names"], results):
            motif_length = len(matrix[0])
            for match in result:
                start = region.start + match.pos
                end = start + motif_length
                score = round(match.score, 5)
                site = OneRegion([region.chrom, start, end, name, score, strand])
                sites.append(site)
        return sites

    # --------------- clustering (optional gimmemotifs) --------------------
    def cluster(self, threshold: float = 0.5, metric: str = "pcc", clust_method: str = "average") -> Dict[str, "MotifList"]:
        """Cluster motifs using gimmemotifs. Returns dict of {cluster_name: MotifList}."""
        try:
            from gimmemotifs.motif import Motif  # noqa: F401
            from gimmemotifs.comparison import MotifComparer
        except Exception as e:  # pragma: no cover
            raise ImportError(
                "Clustering requires gimmemotifs. Install with: `poetry add gimmemotifs`"
            ) from e
        if sns is not None:
            sns.set_style("ticks")
        motif_list = [motif.get_gimmemotif().gimme_obj for motif in self]
        mc = MotifComparer()
        score_dict = mc.get_all_scores(motif_list, motif_list, match="total", metric=metric, combine="mean")
        self.similarity_matrix = generate_similarity_matrix(score_dict)
        import scipy.spatial.distance as ssd
        from scipy.cluster.hierarchy import linkage, fcluster
        vector = ssd.squareform(self.similarity_matrix, checks=not is_symmetric(self.similarity_matrix))
        self.linkage_mat = linkage(vector, method=clust_method)
        fclust_labels = fcluster(self.linkage_mat, threshold, criterion="distance")
        formatted_labels = [f"Cluster_{label}" for label in fclust_labels]
        cluster_dict: Dict[str, MotifList] = {label: MotifList() for label in formatted_labels}
        for i, label in enumerate(formatted_labels):
            cluster_dict[label].append(self[i])
        return cluster_dict

    def create_consensus(self, metric: str = "pcc") -> "OneMotif":
        """Create a consensus OneMotif from this MotifList using gimmemotifs."""
        try:
            from gimmemotifs.motif import Motif
            from gimmemotifs.comparison import MotifComparer
        except Exception as e:  # pragma: no cover
            raise ImportError(
                "Consensus creation requires gimmemotifs. Install with: `poetry add gimmemotifs`"
            ) from e
        self_local = [m.get_gimmemotif() if m.gimme_obj is None else m for m in self]
        motif_list = [m.gimme_obj for m in self_local]
        if len(motif_list) > 1:
            consensus_found = False
            mc = MotifComparer()
            score_dict = mc.get_all_scores(motif_list, motif_list, match="total", metric=metric, combine="mean")
            while not consensus_found:
                best_similarity_motifs = sorted(find_best_pair(motif_list, score_dict))
                new_motif = merge_motifs(motif_list[best_similarity_motifs[0]], motif_list[best_similarity_motifs[1]], metric=metric)
                del motif_list[best_similarity_motifs[1]]
                motif_list[best_similarity_motifs[0]] = new_motif
                if len(motif_list) == 1:
                    consensus_found = True
                else:
                    score_dict[new_motif.id] = score_dict.get(new_motif.id, {})
                    for m in motif_list:
                        score_dict[new_motif.id][m.id] = mc.compare_motifs(new_motif, m, metric=metric)
                        score_dict[m.id][new_motif.id] = mc.compare_motifs(m, new_motif, metric=metric)
        gimmemotif_consensus = motif_list[0]
        gimme_id = gimmemotif_consensus.id
        pwm = [[round(f, 5) for f in l] for l in gimmemotif_consensus.pwm]
        gimmemotif_consensus = Motif(pwm)
        gimmemotif_consensus.id = gimme_id
        onemotif_consensus = gimmemotif_to_onemotif(gimmemotif_consensus)
        onemotif_consensus.gimme_obj = gimmemotif_consensus
        all_names = [m.name for m in self]
        onemotif_consensus.name = ",".join(all_names[:3]) + ("(...)" if len(all_names) > 3 else "")
        return onemotif_consensus

    # --------------- misc --------------------------------------------------
    def plot_motifs(self, nrow=None, ncol=None, output: str = "motif_plot.png", figsize=None, formation: str | List[Tuple[int,int]] = "row"):
        n_motifs = len(self)
        formation, nrow, ncol = get_formation(formation, ncol, nrow, n_motifs)
        if nrow * ncol < n_motifs:
            sys.exit(
                f"ERROR: Insufficient space in grid. motifs={n_motifs}, slots={ncol*nrow}. Increase rows/cols."
            )
        longest_motif = max([len(i[0]) for i in [m.counts for m in self]])
        if figsize is None:
            figsize = (longest_motif * 0.55 * ncol, nrow * 3)
        plt.subplots(squeeze=False, figsize=figsize)
        for x, motif in enumerate(self):
            ax = plt.subplot2grid((nrow, ncol), formation[x])
            motif.create_logo(ax, longest_motif)
        plt.savefig(output)

    def make_unique(self) -> "MotifList":
        seen = {}
        for motif in self:
            m_id = motif.id
            if m_id not in seen:
                seen[m_id] = 1
            else:
                motif.id = f"{motif.id}_{seen[m_id]}"
                seen[m_id] += 1
        return self

    def get_reverse(self) -> "MotifList":
        return MotifList([motif.get_reverse() for motif in self])

    def __add__(self, other):  # pragma: no cover
        return MotifList(list(self) + list(other))


# ----------------------------------------------------------------------------
# conversions & helpers

def gimmemotif_to_onemotif(gimmemotif_obj):
    counts = np.array(gimmemotif_obj.pfm).T.tolist()
    onemotif_obj = OneMotif(motifid=gimmemotif_obj.id, counts=counts)
    return onemotif_obj


def generate_similarity_matrix(score_dict) -> pd.DataFrame:
    m1_keys = list(score_dict.keys())
    m2_keys = list(score_dict.values())[0].keys()
    m1_labels = [s.replace("\t", " ") for s in m1_keys]
    m2_labels = [s.replace("\t", " ") for s in m2_keys]
    similarity_dict = {m: {} for m in m1_labels}
    for i, m1 in enumerate(m1_keys):
        for j, m2 in enumerate(m2_keys):
            score = round(1 - np.mean([score_dict[m1][m2][0], score_dict[m2][m1][0]]), 3)
            score = min(score, 1)
            similarity_dict[m1_labels[i]][m2_labels[j]] = score
            similarity_dict[m2_labels[j]][m1_labels[i]] = score
    similarity_dict_format = {m1: [similarity_dict[m1][m2] for m2 in m2_labels] for m1 in m1_labels}
    df = pd.DataFrame(similarity_dict_format, index=m2_labels).replace(-0, 0)
    return df


def merge_motifs(motif_1, motif_2, metric: str = "pcc"):
    from gimmemotifs.comparison import MotifComparer
    mc = MotifComparer()
    _, pos, orientation = mc.compare_motifs(motif_1, motif_2, metric=metric)
    consensus = motif_1.average_motifs(motif_2, pos=pos, orientation=orientation)
    consensus.id = motif_1.id + "+" + motif_2.id
    return consensus


def find_best_pair(cluster_motifs, score_dict):
    best_similarity = 0
    best_similarity_motifs = [0, 1]
    for i, m in enumerate(cluster_motifs):
        for j, n in enumerate(cluster_motifs):
            if m.id is not n.id:
                this_similarity = score_dict[m.id][n.id][0]
                if this_similarity > best_similarity:
                    best_similarity = this_similarity
                    best_similarity_motifs = [i, j]
    return best_similarity_motifs


def get_formation(formation, ncol, nrow, nmotifs):
    if formation != "alltoone":
        if ncol is None and nrow is None:
            half = int(math.ceil(math.sqrt(nmotifs)))
            ncol, nrow = half, half
        else:
            if ncol is None:
                ncol = int(math.ceil(nmotifs / nrow))
            if nrow is None:
                nrow = int(math.ceil(nmotifs / ncol))
    if isinstance(formation, str):
        if formation == "row":
            formation = []
            rows = list(range(nrow))
            for row in rows:
                for col in range(ncol):
                    formation.append((row, col))
        elif formation == "col":
            formation = []
            rows = list(range(nrow))
            for col in range(ncol):
                for row in rows:
                    formation.append((row, col))
        elif formation == "alltoone":
            formation = []
            rows = list(range(nmotifs - 1))
            for row in rows:
                formation.append((row, 0))
            formation.append((int(math.ceil(len(rows) / 2.0)) - 1, 1))
            ncol = 2
            nrow = len(rows)
        else:
            sys.exit("ERROR: Unknown formation setting.")
    else:
        formation_max_row = max([i[0] for i in formation])
        formation_max_col = max([i[1] for i in formation])
        if nrow < formation_max_row or ncol < formation_max_col:
            sys.exit("ERROR: Grid is too small for specified formation")
    return formation, nrow, ncol


# --------------------------------------------------------------------------
# The core motif object

class OneMotif:
    id: str = ""
    name: str = ""
    bases: List[str] = ["A", "C", "G", "T"]
    bg: np.ndarray = np.array([0.25, 0.25, 0.25, 0.25])
    strands: str = "+ -"
    n: int = 20
    length: Optional[int] = None

    threshold: Optional[float] = None
    gimme_obj = None
    info: Dict[str, str] = {}
    prefix: str = ""

    counts: Optional[List[List[float]]] = None
    pfm: Optional[np.ndarray] = None
    pssm: Optional[np.ndarray] = None

    def __init__(self, motifid: Optional[str], counts: Optional[List[List[float]]] = None, name: Optional[str] = None):
        self.id = motifid if motifid is not None else ""
        self.name = name if name is not None else ""
        if counts is not None:
            self.set_counts(counts)

    def __str__(self):  # pragma: no cover
        return f"{self.__dict__}"

    # ------- naming --------------------------------------------------------
    def set_prefix(self, naming: str = "name_id") -> "OneMotif":
        if naming == "name":
            prefix = self.name
        elif naming == "id":
            prefix = self.id
        elif naming == "name_id":
            prefix = f"{self.name}_{self.id}"
        elif naming == "id_name":
            prefix = f"{self.id}_{self.name}"
        else:
            prefix = "None"
        self.prefix = filafy(prefix)
        return self

    # ------- matrices ------------------------------------------------------
    def get_pfm(self) -> "OneMotif":
        self.pfm = self.counts / np.sum(self.counts, axis=0)  # type: ignore[operator]
        return self

    def get_gimmemotif(self) -> "OneMotif":
        from gimmemotifs.motif import Motif
        self.length = len(self.counts[0])  # type: ignore[index]
        motif_rows = []
        for pos_id in range(self.length):  # type: ignore[arg-type]
            row = [self.counts[letter][pos_id] for letter in range(len(self.bases))]  # type: ignore[index]
            motif_rows.append(row)
        self.gimme_obj = Motif(pfm=motif_rows)
        self.gimme_obj.id = self.id + " " + self.name
        return self

    def get_reverse(self) -> "OneMotif":
        rev_counts = [[], [], [], []]
        rev_counts[0] = self.counts[3][::-1]  # type: ignore[index]
        rev_counts[1] = self.counts[2][::-1]  # type: ignore[index]
        rev_counts[2] = self.counts[1][::-1]  # type: ignore[index]
        rev_counts[3] = self.counts[0][::-1]  # type: ignore[index]
        rev_bg = self.bg[[3, 2, 1, 0]]
        reverse_motif = OneMotif(motifid=self.id, counts=rev_counts, name=self.name)
        reverse_motif.info = self.info
        reverse_motif.bg = rev_bg
        reverse_motif.prefix = self.prefix
        reverse_motif.threshold = self.threshold
        return reverse_motif

    def get_pssm(self, ps: float = 0.01) -> "OneMotif":
        if self.pfm is None:
            self.get_pfm()
        bg_col = self.bg.reshape((-1, 1))
        pseudo_vector = ps * bg_col
        pfm_pc = np.true_divide(self.pfm + pseudo_vector, np.sum(self.pfm + pseudo_vector, axis=0))  # type: ignore[operator]
        self.pssm = np.log(pfm_pc) - np.log(bg_col)
        self.pssm = np.nan_to_num(self.pssm)
        return self

    def get_threshold(self, pvalue: float = 1e-4) -> "OneMotif":
        if self.pssm is None:
            self.get_pssm()
        pssm_tuple = tuple([tuple(row) for row in self.pssm])
        self.threshold = MOODS.tools.threshold_from_p(pssm_tuple, self.bg, pvalue, 4)
        return self

    # ------- metrics -------------------------------------------------------
    def information_content(self, ps: float = 0.01) -> "OneMotif":
        if self.pfm is None:
            self.get_pfm()
        bg_col = self.bg.reshape((-1, 1))
        self.information = self.pfm * (np.log2(self.pfm + ps) - np.log2(bg_col + ps))  # type: ignore[operator]
        self.ic = np.sum(self.information)
        self.bits = self.pfm * np.sum(self.information, axis=0)  # type: ignore[operator]
        return self

    def gc_content(self) -> "OneMotif":
        if self.pfm is None:
            self.get_pfm()
        self.gc_positions = self.pfm[self.bases.index("G")] + self.pfm[self.bases.index("C")]  # type: ignore[index]
        self.gc = np.mean(self.gc_positions)
        return self

    # ------- rendering -----------------------------------------------------
    def logo_to_file(self, filename: str, ylim: Tuple[float, float] = (0, 2)) -> None:
        ext = os.path.splitext(filename)[-1].lower()
        if ext == ".jpg":
            filename = filename[:-4] + ".png"
        logo = self.create_logo(ylim=ylim)
        logo.fig.savefig(filename)  # type: ignore[attr-defined]
        plt.close(logo.fig)

    def get_base(self) -> "OneMotif":
        image = io.BytesIO()
        logo = self.create_logo()
        logo.fig.savefig(image)  # type: ignore[attr-defined]
        self.base = base64.encodebytes(image.getvalue()).decode("utf-8").replace("\n", "")
        return self

    def create_logo(self, ax=None, motif_len: Optional[int] = None, ylim: Tuple[float, float] | str = (0, 2)):
        df = pd.DataFrame(self.counts).transpose()  # type: ignore[arg-type]
        df.columns = self.bases
        if not motif_len:
            motif_len = df.shape[0]
        info_df = logomaker.transform_matrix(df, from_type="counts", to_type="information")
        self.info_df = info_df
        if ylim == "auto":
            ylim = (0, info_df.sum(axis=1).max())
        if not isinstance(ylim, (list, tuple)):
            raise ValueError("ylim should be 'auto' or a tuple of (ymin, ymax)")
        tick_num = 4
        step = (ylim[1] - ylim[0]) / tick_num
        yticks = [ylim[0] + step * i for i in range(tick_num + 1)]
        logo = logomaker.Logo(info_df, ax=ax)
        logo.style_xticks(rotation=0, fmt='%d', anchor=0)
        logo.ax.set_ylim(*ylim)
        logo.ax.set_xlim(-0.5, motif_len - 0.5)
        logo.ax.set_yticks(yticks, minor=False)
        logo.ax.set_xticklabels(range(1, motif_len + 1))
        logo.style_spines(visible=False)
        logo.style_spines(spines=['left', 'bottom'], visible=True)
        logo.ax.xaxis.set_ticks_position('none')
        logo.ax.xaxis.set_tick_params(pad=-1)
        return logo

    # ------- I/O -----------------------------------------------------------
    def set_counts(self, counts: List[List[float]]) -> "OneMotif":
        if len(counts) != 4:
            raise ValueError(
                "Input counts must be of length 4 (ACGT). Got length {0} for motif id='{1}', name='{2}'.".format(
                    len(counts), self.id, self.name
                )
            )
        lengths = [len(base) for base in counts]
        if len(set(lengths)) != 1:
            raise ValueError("All lists in counts must be of the same length.")
        self.counts = counts
        self.length = lengths[0]
        self.n = int(np.sum([row[0] for row in counts]))
        return self

    def as_string(self, output_format: str = "pfm", header: bool = True) -> str:
        out = ""
        if output_format in ["pfm", "jaspar"]:
            out += f">{self.id}\t{self.name}\n"
            for index, base in enumerate(self.bases):
                row = " ".join(map(str, self.counts[index]))
                if output_format == "jaspar":
                    row = f"{base} [{row} ]"
                out += row + "\n"
        elif output_format == "meme":
            if header:
                meme_header = "MEME version 4\n\n"
                meme_header += f"ALPHABET= {''.join(self.bases)}\n\n"
                meme_header += f"strands: {self.strands}\n\n"
                meme_header += "Background letter frequencies\n"
                meme_header += " ".join([f"{self.bases[i]} {self.bg[i]}" for i in range(4)]) + "\n"
                out += meme_header
            out += f"\nMOTIF {self.id} {self.name}\n"
            out += f"letter-probability matrix: alength= {len(self.bases)} w= {self.length} nsites= {int(round(self.n))} E= 0\n"
            if self.pfm is None:
                self.get_pfm()
            precision = 6
            for row in self.pfm.T:
                out += " {0}\n".format("  ".join(map(lambda f: format(round(f, precision), f".{precision}f"), row)))
        else:
            raise ValueError(f"Unsupported motif format: {output_format}")
        return out

    def to_file(self, output_file: str, fmt: str = "pfm") -> "OneMotif":
        with open(output_file, "w") as f:
            f.write(self.as_string(output_format=fmt))
        return self

    @staticmethod
    def from_fasta(fasta: str, motifid: str, name: Optional[str] = None) -> "OneMotif":
        with open(fasta) as handle:
            motif = bio_motifs.read(handle, "sites")
        return OneMotif(
            motifid=motifid,
            counts=[motif.counts[base] for base in ["A", "C", "G", "T"]],
            name=name,
        )

    def __repr__(self) -> str:  # pragma: no cover
        return f"<OneMotif: {self.id}{(' ' + self.name) if len(self.name) > 0 else ''}>"


# --------------------------------------------------------------------------
# format detection

def get_motif_format(content: str) -> str:
    if re.match(r".*MEME version.+", content, re.DOTALL) is not None:
        return "meme"
    if re.match(r">.+A.+\[", content, re.DOTALL) is not None:
        return "jaspar"
    if re.match(r">.+", content, re.DOTALL) is not None:
        return "pfm"
    if re.match(r"AC\s.+", content, re.DOTALL) is not None:
        return "transfac"
    return "unknown"
