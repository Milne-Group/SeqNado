#!/usr/bin/env python3
"""SeqNado CLI application.

This module serves as the entry point for the SeqNado CLI.
Commands are defined in the commands/ submodule and self-register
via @app.command() decorators.
"""

from __future__ import annotations

# Import the app instance (created in app_instance.py)
from seqnado.cli.app_instance import app

# Import command modules to trigger @app.command() registration
# These imports cause each command to register itself with the app
from seqnado.cli.commands import init  # noqa: F401
from seqnado.cli.commands import config  # noqa: F401
from seqnado.cli.commands import tools  # noqa: F401
from seqnado.cli.commands import download  # noqa: F401
from seqnado.cli.commands import design  # noqa: F401
from seqnado.cli.commands import pipeline  # noqa: F401

# Import and register the genomes sub-application
from seqnado.cli.commands.genomes import app as genomes_app

# Add genomes as a sub-application
app.add_typer(genomes_app, name="genomes")

# Optional: prettier tracebacks/console with rich if available
try:
    from rich.traceback import install as _rich_tb_install

    _rich_tb_install(show_locals=False)
except Exception:
    pass


# -------------------------------- Entrypoint --------------------------------
if __name__ == "__main__":
    from seqnado.cli.utils import _resolve_working_dir
    
    _resolve_working_dir()
    app()
