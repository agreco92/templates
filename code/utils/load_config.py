"""
Python equivalent of load_config() from functions.R.

Loads YAML configuration with glue-style interpolation ({projectDir}, {sample_id})
and sample-level overrides, matching the R templates system.

Usage:
    from utils.load_config import load_config

    config = load_config(
        project_dir="/home/user/projects/myproject",
        module_name="02_analysis",
        script_name="04_sciRED_factor_analysis",
        sample_id="sample1"
    )
"""

import os
import re
import copy
import yaml


def _glue_recursive(obj, env):
    """Recursively interpolate {key} placeholders in strings using env dict."""
    if isinstance(obj, dict):
        return {k: _glue_recursive(v, env) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [_glue_recursive(item, env) for item in obj]
    elif isinstance(obj, str):
        try:
            # Replace {key} patterns with values from env
            def replacer(match):
                key = match.group(1)
                if key in env:
                    return str(env[key])
                return match.group(0)  # leave unresolved
            return re.sub(r"\{(\w+)\}", replacer, obj)
        except Exception:
            return obj
    else:
        return obj


def load_config(project_dir, module_name, script_name, sample_id=None):
    """
    Load YAML configuration matching the R load_config() pattern.

    Parameters
    ----------
    project_dir : str
        Project root directory (equivalent to --dir).
    module_name : str
        YAML config file name without .yaml extension.
    script_name : str
        Top-level key in the YAML for this script's defaults.
    sample_id : str, optional
        Sample ID for sample-level overrides.

    Returns
    -------
    dict
        Merged configuration with glue interpolation applied.
    """
    config = {
        "projectDir": project_dir,
        "sample_id": sample_id,
    }

    yaml_path = os.path.join(project_dir, "config", f"{module_name}.yaml")

    if not os.path.exists(yaml_path):
        print(f"WARNING: Config file not found: {yaml_path}")
        return config

    with open(yaml_path) as f:
        yaml_config = yaml.safe_load(f)

    # Load script-level defaults
    if script_name in yaml_config:
        script_config = _glue_recursive(yaml_config[script_name], config)
        config.update(script_config)

    # Load sample-level overrides
    if sample_id and sample_id in yaml_config:
        sample_block = yaml_config[sample_id]
        if script_name in sample_block:
            sample_config = _glue_recursive(sample_block[script_name], config)
            config.update(sample_config)

    return config
