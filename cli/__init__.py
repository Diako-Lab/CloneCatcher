"""
ClusterCatcher CLI Package
==========================

Command-line interface for the ClusterCatcher single-cell analysis pipeline.

Commands:
- sample-information: Process sample CSV and create sample dictionary pickle
- create-config: Generate master config YAML for pipeline  
- run-config: Execute the Snakemake pipeline
"""

__version__ = '1.4.2'

from cli.inspect_libraries import inspect_libraries
from cli.sample_information import sample_information
from cli.run_config import run_config

__all__ = ['inspect_libraries', 'sample_information', 'create_config_cmd', 'run_config']
