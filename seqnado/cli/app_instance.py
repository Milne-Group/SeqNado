"""SeqNado CLI application instance.

This module defines the main Typer app instance and is imported by command modules
to avoid circular dependencies.
"""

from __future__ import annotations

import typer
from loguru import logger

from seqnado._version import __version__
from seqnado.cli.utils import _configure_logging


def version_callback(value: bool):
    """Print version and exit."""
    if value:
        typer.echo(f"SeqNado version {__version__}")
        raise typer.Exit()


# Create the main Typer application
app = typer.Typer(
    add_completion=True,
    no_args_is_help=True,
    help="""
[bold]SeqNado CLI[/bold]

Initialize your environment, build configs, create design files, and run pipelines.
Use --help on any subcommand for details.
""",
)


@app.callback()
def main(
    version: bool = typer.Option(
        False,
        "--version",
        "-v",
        help="Show version and exit.",
        callback=version_callback,
        is_eager=True,
    ),
    verbose: bool = typer.Option(
        False,
        "--verbose",
        help="Enable verbose logging (DEBUG level).",
    ),
):
    """SeqNado CLI main entry point with global options."""
    # Configure logging for all commands
    _configure_logging(verbose)
