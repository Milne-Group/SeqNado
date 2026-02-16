"""List and explore bioinformatics tools command."""
from __future__ import annotations

from typing import Optional

import typer
from loguru import logger

from seqnado.cli.app_instance import app
from seqnado.cli.tools_helpers import _OptionalCategoryCommand
from seqnado.cli.utils import _configure_logging, _RICH_CONSOLE, verbose_option


@app.command(
    cls=_OptionalCategoryCommand,
    help="List and explore bioinformatics tools available in the SeqNado pipeline.",
)
def tools(
    tool_name: Optional[str] = typer.Argument(
        None,
        metavar="[TOOL]",
        help="Specific tool name to get help for (e.g., fastqc, deeptools)",
    ),
    list_tools: bool = typer.Option(
        False, "--list", "-l", help="List all available tools with descriptions."
    ),
    category: Optional[str] = typer.Option(
        None,
        "--category",
        "-c",
        help=(
            "Filter tools by category name or number. Use without a value to "
            "interactively select a category."
        ),
    ),
    show_options: bool = typer.Option(
        False,
        "--options",
        help="Show tool options/help from the container (requires tool name and apptainer).",
    ),
    show_citation: bool = typer.Option(
        False,
        "--citation",
        help="Show the BibTeX citation for a tool (requires tool name).",
    ),
    subcommand: Optional[str] = typer.Option(
        None,
        "--subcommand",
        "-s",
        help="Specify a tool subcommand for help (e.g. plotHeatmap).",
    ),
    verbose: bool = verbose_option(),
) -> None:
    """
    Explore bioinformatics tools in the SeqNado pipeline.

    Args:
        tool_name: Specific tool name to get help for (e.g., fastqc, deeptools)
        list_tools: List all available tools with descriptions
        category: Filter tools by category name or number
        show_options: Show tool options/help from the container (requires tool name and apptainer)
        show_citation: Show the BibTeX citation for a tool (requires tool name)
        subcommand: Specify a tool subcommand for help (e.g. plotHeatmap)
        verbose: Increase logging verbosity

    Examples:
        seqnado tools                         # List all tools
        seqnado tools fastqc                  # Show info for fastqc
        seqnado tools fastqc --options        # Show fastqc --help from container
        seqnado tools fastqc --citation       # Show BibTeX citation for fastqc
        seqnado tools deeptools --options -s plotHeatmap
        seqnado tools --category 5            # List tools in category #5
        seqnado tools -c "Quality Control"    # List tools in a category by name
    """
    from seqnado.tools import (
        format_citation,
        get_tool_citation,
        get_tool_help,
        get_tool_info,
        get_tool_subcommands,
        get_tool_version,
        is_apptainer_available,
        run_tool_help_in_container,
    )
    from seqnado.tools import (
        list_tools as list_all_tools,
    )
    from seqnado.cli.tools_helpers import (
        _validate_subcommand,
        _display_tools_list,
        _resolve_category,
    )

    _configure_logging(verbose)

    # If --options is requested, show tool help from container
    if show_options:
        if not tool_name:
            logger.error("Tool name required with --options")
            raise typer.Exit(code=1)

        if not is_apptainer_available():
            logger.error("Apptainer/Singularity is required to use --options flag")
            logger.info("Please install apptainer or singularity and try again")
            raise typer.Exit(code=1)

        tool_info = get_tool_info(tool_name)
        if not tool_info:
            logger.error(f"Tool '{tool_name}' not found")
            logger.info("Run 'seqnado tools --list' to see all available tools.")
            raise typer.Exit(code=1)

        if subcommand:
            subcommand = _validate_subcommand(
                tool_name, subcommand, get_tool_subcommands(tool_name)
            )

        logger.info(f"Fetching help for '{tool_name}' from SeqNado container...")
        help_output = run_tool_help_in_container(tool_name, subcommand=subcommand)
        if _RICH_CONSOLE:
            _RICH_CONSOLE.print(help_output)
        else:
            typer.echo(help_output)
        raise typer.Exit(code=0)

    # If --citation is requested, show BibTeX entry
    if show_citation:
        if not tool_name:
            logger.error("Tool name required with --citation")
            raise typer.Exit(code=1)

        tool_info = get_tool_info(tool_name)
        if not tool_info:
            logger.error(f"Tool '{tool_name}' not found")
            logger.info("Run 'seqnado tools --list' to see all available tools.")
            raise typer.Exit(code=1)

        bibtex = get_tool_citation(tool_name)
        if bibtex:
            formatted = format_citation(bibtex)
            if _RICH_CONSOLE:
                if formatted:
                    _RICH_CONSOLE.print(f"\n[bold cyan]Citation for {tool_name}:[/bold cyan]\n")
                    _RICH_CONSOLE.print(formatted)
                _RICH_CONSOLE.print(f"\n[bold cyan]BibTeX:[/bold cyan]\n")
                _RICH_CONSOLE.print(bibtex)
            else:
                if formatted:
                    typer.echo(f"\nCitation for {tool_name}:\n")
                    typer.echo(formatted)
                typer.echo("\nBibTeX:\n")
                typer.echo(bibtex)
        else:
            logger.warning(f"No citation available for '{tool_name}'")
        raise typer.Exit(code=0)

    # Resolve category (number -> name)
    if category:
        category = _resolve_category(category)

    # List tools (default when no tool specified, or with --list / --category)
    if list_tools or tool_name is None:
        tools_list = list_all_tools(category=category)
        _display_tools_list(tools_list, category=category)
        raise typer.Exit(code=0)

    # Show help for specific tool
    tool_info = get_tool_info(tool_name)

    if not tool_info:
        logger.error(f"Tool '{tool_name}' not found.")
        logger.info("Run 'seqnado tools --list' to see all available tools.")
        raise typer.Exit(code=1)

    subcommands = get_tool_subcommands(tool_name)
    if subcommand:
        subcommand = _validate_subcommand(tool_name, subcommand, subcommands)

    # Display tool information
    if _RICH_CONSOLE:
        _RICH_CONSOLE.print(f"\n[bold cyan]Tool: {tool_name}[/bold cyan]")
        _RICH_CONSOLE.print(f"[bold]Description:[/bold] {tool_info['description']}")
        _RICH_CONSOLE.print(f"[bold]Category:[/bold] {tool_info['category']}")
        _RICH_CONSOLE.print(f"[bold]Command:[/bold] {tool_info['command']}")
        if subcommands:
            _RICH_CONSOLE.print(
                f"[bold]Subcommands:[/bold] {', '.join(subcommands)}"
            )

        # Try to get version info
        _RICH_CONSOLE.print("\n[bold cyan]Version Information:[/bold cyan]")
        version_info = get_tool_version(
            tool_name,
            use_container=is_apptainer_available(),
            subcommand=subcommand,
        )
        _RICH_CONSOLE.print(f"[bold white]Version: {version_info}[/bold white]")

        # Get and display help
        _RICH_CONSOLE.print("\n[bold cyan]Help / Usage:[/bold cyan]")
        help_text = get_tool_help(tool_name, subcommand=subcommand)
        if (
            help_text.startswith("Help information not available")
            and is_apptainer_available()
        ):
            help_text = run_tool_help_in_container(tool_name, subcommand=subcommand)
        _RICH_CONSOLE.print(help_text)

        if is_apptainer_available():
            _RICH_CONSOLE.print(
                "\n[dim]Tip: Use 'seqnado tools "
                + tool_name
                + " --options' to see detailed options from the container.[/dim]"
            )
    else:
        typer.echo(f"\nTool: {tool_name}")
        typer.echo(f"Description: {tool_info['description']}")
        typer.echo(f"Category: {tool_info['category']}")
        typer.echo(f"Command: {tool_info['command']}")
        if subcommands:
            typer.echo(f"Subcommands: {', '.join(subcommands)}")
        typer.echo("\nVersion Information:")
        version_info = get_tool_version(
            tool_name,
            use_container=is_apptainer_available(),
            subcommand=subcommand,
        )
        typer.echo(version_info)
        typer.echo("\nHelp / Usage:")
        help_text = get_tool_help(tool_name, subcommand=subcommand)
        if (
            help_text.startswith("Help information not available")
            and is_apptainer_available()
        ):
            help_text = run_tool_help_in_container(tool_name, subcommand=subcommand)
        typer.echo(help_text)
