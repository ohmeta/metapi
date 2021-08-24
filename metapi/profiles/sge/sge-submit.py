#!/usr/bin/env python3

import os
import re
import math
import argparse
import subprocess

# use warnings.warn() rather than print() to output info in this script
# because snakemake expects the jobid to be the only output
import warnings

from snakemake import io
from snakemake.utils import read_job_properties

DEFAULT_JOB_NAME = "snakemake_job"
QSUB_DEFAULTS = "-cwd -V"
CLUSTER_CONFIG = "cluster.yaml"

# SGE syntax for options is `-option [value]` and for resources is `-l name=value`
# we therefore distinguish the two in this script to make it easier to handle.
# We also define some aliases for options and resources so that the rules can
# be more expressive than a list of cryptic SGE resources.

# We additionally pickup a list of environment modules which will be loaded in the
# jobscript

OPTION_MAPPING = {
    "binding": ("binding",),
    "cwd"    : ("cwd",),
    "e"      : ("e", "error"),
    "hard"   : ("hard",),
    "j"      : ("j", "join"),
    "m"      : ("m", "mail_options"),
    "M"      : ("M", "email"),
    "notify" : ("notify",),
    "now"    : ("now",),
    "N"      : ("N", "name"),
    "o"      : ("o", "output"),
    "P"      : ("P", "project"),
    "p"      : ("p", "priority"),
    "pe"     : ("pe", "parallel_environment"),
    "pty"    : ("pty",),
    "q"      : ("q", "queue"),
    "R"      : ("R", "reservation"),
    "r"      : ("r", "rerun"),
    "soft"   : ("soft",),
    "v"      : ("v", "variable"),
    "V"      : ("V", "export_env")
}

RESOURCE_MAPPING = {
    # default queue resources
    "qname"            : ("qname",),
    "hostname"         : ("hostname",),
    # "notify" -- conflicts with OPTION_MAPPING
    "calendar"         : ("calendar",),
    "min_cpu_interval" : ("min_cpu_interval",),
    "tmpdir"           : ("tmpdir",),
    "seq_no"           : ("seq_no",),
    "s_rt"             : ("s_rt", "soft_runtime", "soft_walltime"),
    "h_rt"             : ("h_rt", "time", "runtime", "walltime"),
    "s_cpu"            : ("s_cpu", "soft_cpu"),
    "h_cpu"            : ("h_cpu", "cpu"),
    "s_data"           : ("s_data", "soft_data"),
    "h_data"           : ("h_data", "data"),
    "s_stack"          : ("s_stack", "soft_stack"),
    "h_stack"          : ("h_stack", "stack"),           
    "s_core"           : ("s_core", "soft_core"),
    "h_core"           : ("h_core", "core"),
    "s_rss"            : ("s_rss", "soft_resident_set_size"),
    "h_rss"            : ("h_rss", "resident_set_size"),
    # default host resources
    "slots"            : ("slots",),
    "s_vmem"           : ("s_vmem", "soft_memory", "soft_virtual_memory"),
    "h_vmem"           : ("h_vmem", "mem", "memory", "virtual_memory"),
    "s_fsize"          : ("s_fsize", "soft_file_size"),
    "h_fsize"          : ("h_fsize", "file_size"),
}

def add_custom_resources(resources, resource_mapping=RESOURCE_MAPPING):
    """Adds new resources to resource_mapping.

       resources -> dict where key is sge resource name and value is a 
                    single name or a list of names to be used as aliased
    """
    for key, val in resources.items():
        if key not in resource_mapping:
            resource_mapping[key] = tuple()

        # make sure the resource name itself is an alias
        resource_mapping[key] += (key,)
        if isinstance(val, list):
            for alias in val:
                if val != key:
                    resource_mapping[key] += (alias,)
        else:
            if val != key:
                resource_mapping[key] += (val,)

def add_custom_options(options, option_mapping=OPTION_MAPPING):
    """Adds new options to option_mapping.

       options -> dict where key is sge option name and value is a single name
                  or a list of names to be used as aliased
    """
    for key, val in options.items():
        if key not in option_mapping:
            option_mapping[key] = tuple()

        # make sure the option name itself is an alias
        option_mapping[key] += (key,)
        if isinstance(val, list):
            for alias in val:
                if val != key:
                    option_mapping[key] += (alias,)
        else:
            if val != key:
                option_mapping[key] += (val,)

def parse_jobscript():
    """Minimal CLI to require/only accept single positional argument."""
    p = argparse.ArgumentParser(description="SGE snakemake submit script")
    p.add_argument("jobscript", help="Snakemake jobscript with job properties.")
    return p.parse_args().jobscript

def parse_qsub_defaults(parsed):
    """Unpack QSUB_DEFAULTS."""
    d = parsed.split() if type(parsed) == str else parsed
    
    options={}
    for arg in d:
        if "=" in arg:
            k,v = arg.split("=")
            options[k.strip("-")] = v.strip()
        else:
            options[arg.strip("-")] = ""
    return options

def format_job_properties(string):
    # we use 'rulename' rather than 'rule' for consistency with the --cluster-config 
    # snakemake option
    return string.format(rulename=job_properties['rule'], jobid=job_properties['jobid'])


def parse_qsub_settings(source, resource_mapping=RESOURCE_MAPPING, option_mapping=OPTION_MAPPING):
    job_options = { "options" : {}, "resources" : {}}

    for skey, sval in source.items():
        found = False
        for rkey, rval in resource_mapping.items():
            if skey in rval:
                found = True
                # Snakemake resources can only be defined as integers, but SGE interprets
                # plain integers for memory as bytes. This hack means we interpret memory
                # requests as gigabytes
                if (rkey == 's_vmem') or (rkey == 'h_vmem'):
                    job_options["resources"].update({rkey : str(sval) + 'G'})
                else:
                    job_options["resources"].update({rkey : sval})
                break
        if found: continue
        for okey, oval in option_mapping.items():
            if skey in oval:
                found = True
                job_options["options"].update({okey : sval})
                break
        if not found:
            raise KeyError(f"Unknown SGE option or resource: {skey}")

    return job_options

def load_cluster_config(path):
    """Load config to dict either from absolute path or relative to profile dir."""
    if path:
        path = os.path.join(os.path.dirname(__file__), os.path.expandvars(path))
        default_cluster_config = io.load_configfile(path)
    else:
        default_cluster_config = {}
    if "__default__" not in default_cluster_config:
        default_cluster_config["__default__"] = {}
    return default_cluster_config

def ensure_directory_exists(path):
    """Check if directory exists and create if not"""
    directory = os.path.dirname(path)
    if not os.path.exists(directory):
        os.makedirs(directory, exist_ok=True)
    return


def update_double_dict(outer, inner):
    """Similar to dict.update() but does the update on nested dictionaries"""
    for k, v in outer.items():
        outer[k].update(inner[k])

def sge_option_string(key, val):
    if val == "":
        return f"-{key}"
    if type(val) == bool:
        return f"-{key} " + ("yes" if val else "no")
    return format_job_properties(f"-{key} {val}")

def sge_resource_string(key, val):
    if val == "":
        return f"-l {key}"
    if type(val) == bool:
        return f"-{key}=" + ("true" if val else "false")
    return f"-l {key}={val}"

def submit_job(jobscript, qsub_settings):
    """Submit jobscript and return jobid."""
    flatten = lambda l: [item for sublist in l for item in sublist]
    batch_options = flatten([sge_option_string(k,v).split() for k, v in qsub_settings["options"].items()])
    batch_resources = flatten([sge_resource_string(k, v).split() for k, v in qsub_settings["resources"].items()])
    try:
        # -terse means only the jobid is returned rather than the normal 'Your job...' string
        jobid = subprocess.check_output(["qsub", "-terse"] + batch_options + batch_resources + [jobscript]).decode().rstrip()
    except subprocess.CalledProcessError as e:
        raise e
    except Exception as e:
        raise e
    return jobid

qsub_settings = { "options" : {}, "resources" : {}}

jobscript = parse_jobscript()

# get the job properties dictionary from snakemake 
job_properties = read_job_properties(jobscript)

# load the default cluster config
cluster_config = load_cluster_config(CLUSTER_CONFIG)

add_custom_resources(cluster_config["__resources__"])

add_custom_options(cluster_config["__options__"])

# qsub default arguments
update_double_dict(qsub_settings, parse_qsub_settings(parse_qsub_defaults(QSUB_DEFAULTS)))

# cluster_config defaults
update_double_dict(qsub_settings, parse_qsub_settings(cluster_config["__default__"]))

# resources defined in the snakemake file (note that these must be integer)
# we pass an empty dictionary for option_mapping because options should not be
# specified in the snakemake file
update_double_dict(qsub_settings, parse_qsub_settings(job_properties.get("resources", {}), option_mapping={}))

# get any rule specific options/resources from the default cluster config
update_double_dict(qsub_settings, parse_qsub_settings(cluster_config.get(job_properties.get("rule"), {})))

# get any options/resources specified through the --cluster-config command line argument
update_double_dict(qsub_settings, parse_qsub_settings(job_properties.get("cluster", {})))

# ensure qsub output dirs exist
for o in ("o", "e"):
    ensure_directory_exists(qsub_settings["options"][o]) if o in qsub_settings["options"] else None

# submit job and echo id back to Snakemake (must be the only stdout)
print(submit_job(jobscript, qsub_settings))

