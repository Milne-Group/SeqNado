"""
MkDocs hooks for pre-build automation tasks.
"""

import subprocess
import sys
from pathlib import Path


def on_pre_build(config):
    """
    Run before the build starts.
    Regenerates tool citations from the BibTeX file.
    """
    print("Running pre-build hooks...")
    
    # Get the root directory of the project
    docs_dir = Path(config['docs_dir'])
    root_dir = docs_dir.parent
    
    # Run the citation generator script
    script_path = docs_dir / "scripts" / "generate_tool_citations.py"
    
    if script_path.exists():
        print(f"Generating tool citations from BibTeX file...")
        try:
            result = subprocess.run(
                [sys.executable, str(script_path), "--update"],
                cwd=str(root_dir),
                capture_output=True,
                text=True,
                check=True
            )
            print(result.stdout)
            if result.stderr:
                print(result.stderr, file=sys.stderr)
        except subprocess.CalledProcessError as e:
            print(f"Error generating citations: {e}", file=sys.stderr)
            print(e.stdout)
            print(e.stderr, file=sys.stderr)
            # Don't fail the build, just warn
    else:
        print(f"Warning: Citation generator script not found at {script_path}")
    
    print("Pre-build hooks completed.")
