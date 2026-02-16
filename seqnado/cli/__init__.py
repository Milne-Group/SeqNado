"""SeqNado Command Line Interface package.

This package contains the CLI application and related helper modules.

Public API:
    app: The main Typer application (with all commands registered)
"""

from __future__ import annotations

# Import from app.py (not app_instance.py) to trigger command registration
from seqnado.cli.app import app

__all__ = [
    "app",
]
