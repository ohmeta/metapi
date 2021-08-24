#!/usr/bin/env python3
import os
import sys
from os.path import dirname
import re
import math
import argparse
import subprocess as sp
from io import StringIO

from snakemake import io
from snakemake.io import Wildcards
from snakemake.utils import SequenceFormatter
from snakemake.utils import AlwaysQuotedFormatter
from snakemake.utils import QuotedFormatter
from snakemake.exceptions import WorkflowError
from snakemake.logging import logger

from CookieCutter import CookieCutter


def _convert_units_to_mb(memory):
    """If memory is specified with SI unit, convert to MB"""
    if isinstance(memory, int) or isinstance(memory, float):
        return int(memory)
    siunits = {"K": 1e-3, "M": 1, "G": 1e3, "T": 1e6}
    regex = re.compile(r"(\d+)({})$".format("|".join(siunits.keys())))
    m = regex.match(memory)
    if m is None:
        logger.error(
            (
                f"unsupported memory specification '{memory}';"
                "  allowed suffixes: [K|M|G|T]"
            )
        )
        sys.exit(1)
    factor = siunits[m.group(2)]
    return int(int(m.group(1)) * factor)


def parse_jobscript():
    """Minimal CLI to require/only accept single positional argument."""
    p = argparse.ArgumentParser(description="SLURM snakemake submit script")
    p.add_argument("jobscript", help="Snakemake jobscript with job properties.")
    return p.parse_args().jobscript


def parse_sbatch_defaults(parsed):
    """Unpack SBATCH_DEFAULTS."""
    d = parsed.split() if type(parsed) == str else parsed
    args = {}
    for keyval in [a.split("=") for a in d]:
        k = keyval[0].strip().strip("-")
        v = keyval[1].strip() if len(keyval) == 2 else None
        args[k] = v
    return args


def load_cluster_config(path):
    """Load config to dict

    Load configuration to dict either from absolute path or relative
    to profile dir.
    """
    if path:
        path = os.path.join(dirname(__file__), os.path.expandvars(path))
        dcc = io.load_configfile(path)
    else:
        dcc = {}
    if "__default__" not in dcc:
        dcc["__default__"] = {}
    return dcc


# adapted from format function in snakemake.utils
def format(_pattern, _quote_all=False, **kwargs):  # noqa: A001
    """Format a pattern in Snakemake style.
    This means that keywords embedded in braces are replaced by any variable
    values that are available in the current namespace.
    """
    fmt = SequenceFormatter(separator=" ")
    if _quote_all:
        fmt.element_formatter = AlwaysQuotedFormatter()
    else:
        fmt.element_formatter = QuotedFormatter()
    try:
        return fmt.format(_pattern, **kwargs)
    except KeyError as ex:
        raise NameError(
            f"The name {ex} is unknown in this context. Please "
            "make sure that you defined that variable. "
            "Also note that braces not used for variable access "
            "have to be escaped by repeating them "
        )


#  adapted from Job.format_wildcards in snakemake.jobs
def format_wildcards(string, job_properties):
    """ Format a string with variables from the job. """

    class Job(object):
        def __init__(self, job_properties):
            for key in job_properties:
                setattr(self, key, job_properties[key])

    job = Job(job_properties)
    if "params" in job_properties:
        job._format_params = Wildcards(fromdict=job_properties["params"])
    else:
        job._format_params = None
    if "wildcards" in job_properties:
        job._format_wildcards = Wildcards(fromdict=job_properties["wildcards"])
    else:
        job._format_wildcards = None
    _variables = dict()
    _variables.update(
        dict(params=job._format_params, wildcards=job._format_wildcards)
    )
    if hasattr(job, "rule"):
        _variables.update(dict(rule=job.rule))
    try:
        return format(string, **_variables)
    except NameError as ex:
        raise WorkflowError(
            "NameError with group job {}: {}".format(job.jobid, str(ex))
        )
    except IndexError as ex:
        raise WorkflowError(
            "IndexError with group job {}: {}".format(job.jobid, str(ex))
        )


# adapted from ClusterExecutor.cluster_params function in snakemake.executor
def format_values(dictionary, job_properties):
    formatted = dictionary.copy()
    for key, value in list(formatted.items()):
        if key == "mem":
            value = str(_convert_units_to_mb(value))
        if isinstance(value, str):
            try:
                formatted[key] = format_wildcards(value, job_properties)
            except NameError as e:
                msg = "Failed to format cluster config " "entry for job {}.".format(
                    job_properties["rule"]
                )
                raise WorkflowError(msg, e)
    return formatted


def convert_job_properties(job_properties, resource_mapping=None):
    options = {}
    if resource_mapping is None:
        resource_mapping = {}
    resources = job_properties.get("resources", {})
    for k, v in resource_mapping.items():
        options.update({k: resources[i] for i in v if i in resources})

    if "threads" in job_properties:
        options["cpus-per-task"] = job_properties["threads"]
    return options


def ensure_dirs_exist(path):
    """Ensure output folder for Slurm log files exist."""
    di = dirname(path)
    if di == "":
        return
    if not os.path.exists(di):
        os.makedirs(di, exist_ok=True)
    return


def format_sbatch_options(**sbatch_options):
    """Format sbatch options"""
    options = []
    for k, v in sbatch_options.items():
        val = ""
        if v is not None:
            val = f"={v}"
        options.append(f"--{k}{val}")
    return options


def submit_job(jobscript, **sbatch_options):
    """Submit jobscript and return jobid."""
    options = format_sbatch_options(**sbatch_options)
    try:
        cmd = ["sbatch"] + ["--parsable"] + options + [jobscript]
        res = sp.check_output(cmd)
    except sp.CalledProcessError as e:
        raise e
    # Get jobid
    res = res.decode()
    try:
        jobid = re.search(r"(\d+)", res).group(1)
    except Exception as e:
        raise e
    return jobid


def advanced_argument_conversion(arg_dict):
    """Experimental adjustment of sbatch arguments to the given or default partition."""
    # Currently not adjusting for multiple node jobs
    nodes = int(arg_dict.get("nodes", 1))
    if nodes > 1:
        return arg_dict
    partition = arg_dict.get("partition", None) or _get_default_partition()
    constraint = arg_dict.get("constraint", None)
    ncpus = int(arg_dict.get("cpus-per-task", 1))
    runtime = arg_dict.get("time", None)
    memory = _convert_units_to_mb(arg_dict.get("mem", 0))
    config = _get_cluster_configuration(partition, constraint, memory)
    mem = arg_dict.get("mem", ncpus * min(config["MEMORY_PER_CPU"]))
    mem = _convert_units_to_mb(mem)
    if mem > max(config["MEMORY"]):
        logger.info(
            f"requested memory ({mem}) > max memory ({max(config['MEMORY'])}); "
            "adjusting memory settings"
        )
        mem = max(config["MEMORY"])

    # Calculate available memory as defined by the number of requested
    # cpus times memory per cpu
    AVAILABLE_MEM = ncpus * min(config["MEMORY_PER_CPU"])
    # Add additional cpus if memory is larger than AVAILABLE_MEM
    if mem > AVAILABLE_MEM:
        logger.info(
            f"requested memory ({mem}) > "
            f"ncpus x MEMORY_PER_CPU ({AVAILABLE_MEM}); "
            "trying to adjust number of cpus up"
        )
        ncpus = int(math.ceil(mem / min(config["MEMORY_PER_CPU"])))
    if ncpus > max(config["CPUS"]):
        logger.info(
            f"ncpus ({ncpus}) > available cpus ({max(config['CPUS'])}); "
            "adjusting number of cpus down"
        )
        ncpus = min(int(max(config["CPUS"])), ncpus)
    adjusted_args = {"mem": int(mem), "cpus-per-task": ncpus}

    # Update time. If requested time is larger than maximum allowed time, reset
    if runtime:
        runtime = time_to_minutes(runtime)
        time_limit = max(config["TIMELIMIT_MINUTES"])
        if runtime > time_limit:
            logger.info(
                f"time (runtime) > time limit {time_limit}; " "adjusting time down"
            )
            adjusted_args["time"] = time_limit

    # update and return
    arg_dict.update(adjusted_args)
    return arg_dict


timeformats = [
    re.compile(r"^(?P<days>\d+)-(?P<hours>\d+):(?P<minutes>\d+):(?P<seconds>\d+)$"),
    re.compile(r"^(?P<days>\d+)-(?P<hours>\d+):(?P<minutes>\d+)$"),
    re.compile(r"^(?P<days>\d+)-(?P<hours>\d+)$"),
    re.compile(r"^(?P<hours>\d+):(?P<minutes>\d+):(?P<seconds>\d+)$"),
    re.compile(r"^(?P<minutes>\d+):(?P<seconds>\d+)$"),
    re.compile(r"^(?P<minutes>\d+)$"),
]


def time_to_minutes(time):
    """Convert time string to minutes.

    According to slurm:

      Acceptable time formats include "minutes", "minutes:seconds",
      "hours:minutes:seconds", "days-hours", "days-hours:minutes"
      and "days-hours:minutes:seconds".

    """
    if not isinstance(time, str):
        time = str(time)
    d = {"days": 0, "hours": 0, "minutes": 0, "seconds": 0}
    regex = list(filter(lambda regex: regex.match(time) is not None, timeformats))
    if len(regex) == 0:
        return
    assert len(regex) == 1, "multiple time formats match"
    m = regex[0].match(time)
    d.update(m.groupdict())
    minutes = (
        int(d["days"]) * 24 * 60
        + int(d["hours"]) * 60
        + int(d["minutes"])
        + math.ceil(int(d["seconds"]) / 60)
    )
    assert minutes > 0, "minutes has to be greater than 0"
    return minutes


def _get_default_partition():
    """Retrieve default partition for cluster"""
    cluster = CookieCutter.get_cluster_option()
    cmd = f"sinfo -O partition {cluster}"
    res = sp.check_output(cmd.split())
    m = re.search(r"(?P<partition>\S+)\*", res.decode(), re.M)
    partition = m.group("partition")
    return partition


def _get_cluster_configuration(partition, constraints=None, memory=0):
    """Retrieve cluster configuration.

    Retrieve cluster configuration for a partition filtered by
    constraints, memory and cpus

    """
    try:
        import pandas as pd
    except ImportError:
        print(
            "Error: currently advanced argument conversion "
            "depends on 'pandas'.", file=sys.stderr
        )
        sys.exit(1)

    if constraints:
        constraint_set = set(constraints.split(","))
    cluster = CookieCutter.get_cluster_option()
    cmd = f"sinfo -e -o %all -p {partition} {cluster}".split()
    try:
        output = sp.Popen(" ".join(cmd), shell=True, stdout=sp.PIPE).communicate()
    except Exception as e:
        print(e)
        raise
    data = re.sub("^CLUSTER:.+\n", "", re.sub(" \\|", "|", output[0].decode()))
    df = pd.read_csv(StringIO(data), sep="|")
    try:
        df["TIMELIMIT_MINUTES"] = df["TIMELIMIT"].apply(time_to_minutes)
        df["MEMORY_PER_CPU"] = df["MEMORY"] / df["CPUS"]
        df["FEATURE_SET"] = df["AVAIL_FEATURES"].str.split(",").apply(set)
    except Exception as e:
        print(e)
        raise
    if constraints:
        constraint_set = set(constraints.split(","))
        i = df["FEATURE_SET"].apply(lambda x: len(x.intersection(constraint_set)) > 0)
        df = df.loc[i]
    memory = min(_convert_units_to_mb(memory), max(df["MEMORY"]))
    df = df.loc[df["MEMORY"] >= memory]
    return df
