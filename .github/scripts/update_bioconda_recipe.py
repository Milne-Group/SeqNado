#!/usr/bin/env python3
"""
Update bioconda meta.yaml with new version, SHA256, and dependencies.

This script handles the complexity of updating meta.yaml while preserving
Jinja2 templating that bioconda recipes use.

Usage:
    python3 update_bioconda_recipe.py <version> <sha256> <workflow_repo_path>
"""

import sys
import yaml
import re
from pathlib import Path


def extract_jinja2_and_yaml(raw_content):
    """
    Separate Jinja2 statements from YAML content.
    
    Returns:
        tuple: (jinja_lines, yaml_content)
    """
    jinja_lines = []
    yaml_lines = []
    jinja_pattern = r'{%\s*set\s+(\w+)\s*=\s*"([^"]*?)"\s*%}'
    
    for line in raw_content.split('\n'):
        if re.match(jinja_pattern, line.strip()):
            jinja_lines.append(line)
        else:
            yaml_lines.append(line)
    
    return jinja_lines, '\n'.join(yaml_lines)


def replace_jinja_expressions(yaml_content):
    """
    Replace Jinja2 expressions {{ ... }} with placeholders for safe YAML parsing.
    
    Returns:
        tuple: (modified_yaml, jinja_placeholders_dict)
    """
    jinja_placeholders = {}
    placeholder_counter = [0]  # Use list to allow mutation in nested function
    
    def replace_expr(match):
        placeholder = f"__JINJA_{placeholder_counter[0]}__"
        jinja_placeholders[placeholder] = match.group(0)
        placeholder_counter[0] += 1
        return placeholder
    
    modified = re.sub(r'{{.*?}}', replace_expr, yaml_content, flags=re.DOTALL)
    return modified, jinja_placeholders


def restore_jinja_expressions(text, jinja_placeholders):
    """Restore Jinja2 expressions from placeholders."""
    result = text
    for placeholder, original in jinja_placeholders.items():
        result = result.replace(placeholder, original)
    return result


def update_meta_yaml(version, sha256, meta_deps_path):
    """
    Update meta.yaml with new version, SHA256, and dependencies.
    
    Args:
        version: New version string
        sha256: New SHA256 hash
        meta_deps_path: Path to meta_deps.txt containing conda dependencies
    """
    meta_file = Path("meta.yaml")
    
    if not meta_file.exists():
        raise FileNotFoundError(f"{meta_file} not found in {Path.cwd()}")
    
    try:
        meta_deps_exists = Path(meta_deps_path).exists()
    except (OSError, PermissionError):
        meta_deps_exists = False

    if not meta_deps_exists:
        raise FileNotFoundError(f"{meta_deps_path} not found")
    
    # Read conda dependencies
    with open(meta_deps_path) as f:
        conda_deps_raw = f.read()
    
    run_deps = []
    for line in conda_deps_raw.strip().split('\n'):
        dep = line.strip()
        if dep.startswith("- "):
            dep = dep[2:]
        if dep:
            run_deps.append(dep)
    
    # Load meta.yaml and handle Jinja2
    with open(meta_file) as f:
        raw_content = f.read()
    
    jinja_lines, yaml_content = extract_jinja2_and_yaml(raw_content)
    yaml_content, jinja_placeholders = replace_jinja_expressions(yaml_content)
    
    # Parse YAML safely
    meta = yaml.safe_load(yaml_content)
    
    # Update content
    meta["package"]["version"] = version
    meta["source"]["sha256"] = sha256
    meta["requirements"]["run"] = run_deps
    
    # Write back with Jinja2 preserved
    with open(meta_file, "w") as f:
        # Write Jinja2 statements
        for line in jinja_lines:
            line = re.sub(r'version\s*=\s*"[^"]*"', f'version = "{version}"', line)
            f.write(line + '\n')
        
        if jinja_lines:
            f.write('\n')
        
        # Dump YAML and restore Jinja2 expressions
        yaml_str = yaml.dump(meta, default_flow_style=False, sort_keys=False)
        yaml_str = restore_jinja_expressions(yaml_str, jinja_placeholders)
        f.write(yaml_str)
    
    print(f"✓ Updated version to {version}")
    print(f"✓ Updated SHA256")
    print(f"✓ Updated {len(run_deps)} run dependencies")
    print(f"\n=== Updated meta.yaml (first 50 lines) ===")
    
    with open(meta_file) as f:
        for i, line in enumerate(f, 1):
            if i <= 50:
                print(line, end="")
            else:
                break


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: update_bioconda_recipe.py <version> <sha256> <meta_deps_path>")
        sys.exit(1)
    
    version = sys.argv[1]
    sha256 = sys.argv[2]
    meta_deps_path = sys.argv[3]
    
    try:
        update_meta_yaml(version, sha256, meta_deps_path)
    except Exception as e:
        print(f"❌ Error: {e}", file=sys.stderr)
        sys.exit(1)
