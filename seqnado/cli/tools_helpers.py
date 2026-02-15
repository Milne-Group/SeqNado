"""Helper functions and classes for the tools command."""

from __future__ import annotations

import sys
from typing import List, Optional

import click
import typer
from loguru import logger

from seqnado.cli.utils import _RICH_CONSOLE


def _validate_subcommand(
    tool_name: str,
    subcommand: str,
    subcommands: List[str],
) -> str:
    """
    Validate and resolve a subcommand for a tool, exiting on error.
    
    Args:
        tool_name: Name of the tool
        subcommand: Subcommand name to validate
        subcommands: List of valid subcommands for this tool
    
    Returns:
        The validated subcommand name
    
    Raises:
        typer.Exit: If subcommand is invalid
    """
    if not subcommands:
        logger.error(f"Tool '{tool_name}' does not define subcommands.")
        raise typer.Exit(code=1)

    for cmd in subcommands:
        if cmd.lower() == subcommand.lower():
            return cmd

    logger.error(f"Unknown subcommand '{subcommand}' for tool '{tool_name}'.")
    logger.info(f"Available subcommands: {', '.join(subcommands)}")
    raise typer.Exit(code=1)


def _display_tools_list(
    tools_list: list,
    category: Optional[str] = None,
) -> None:
    """
    Display a formatted list of tools grouped by category.
    
    Args:
        tools_list: List of (tool_name, description, category) tuples
        category: Filter by this category (for logging purposes)
    
    Raises:
        typer.Exit: If no tools found
    """
    from seqnado.tools import get_categories

    if _RICH_CONSOLE:
        _RICH_CONSOLE.print(
            "[bold cyan]Available Tools in SeqNado Pipeline[/bold cyan]\n"
        )
    else:
        typer.echo("Available Tools in SeqNado Pipeline\n")

    if not tools_list:
        typer.echo(f"No tools found for category '{category}'")
        raise typer.Exit(code=1)

    tools_by_category: dict = {}
    for tool_name_item, description, cat in tools_list:
        if cat not in tools_by_category:
            tools_by_category[cat] = []
        tools_by_category[cat].append((tool_name_item, description))

    ordered_categories = get_categories()
    for idx, cat in enumerate(ordered_categories, 1):
        if cat not in tools_by_category:
            continue
        if _RICH_CONSOLE:
            _RICH_CONSOLE.print(f"[bold]{idx}. {cat}[/bold]")
            for tool_name_item, description in tools_by_category[cat]:
                _RICH_CONSOLE.print(
                    f"  [cyan]{tool_name_item:<20}[/cyan] {description}"
                )
            _RICH_CONSOLE.print()
        else:
            typer.echo(f"{idx}. {cat}")
            for tool_name_item, description in tools_by_category[cat]:
                typer.echo(f"  {tool_name_item:<20} {description}")
            typer.echo()

    if _RICH_CONSOLE:
        _RICH_CONSOLE.print(
            "Use 'seqnado tools TOOL' to see help for a specific tool."
        )
        _RICH_CONSOLE.print(
            "Use 'seqnado tools TOOL --options' to see tool options from the container."
        )
    else:
        typer.echo("Use 'seqnado tools TOOL' to see help for a specific tool.")
        typer.echo(
            "Use 'seqnado tools TOOL --options' to see tool options from the container."
        )


def _resolve_category(category: str) -> str:
    """
    Resolve a category selection to a category name.

    Handles three forms:
      - ``"__PROMPT__"`` – interactive TTY prompt
      - a digit string – looked up by 1-based index
      - anything else – returned as-is (assumed to be a category name)
    
    Args:
        category: Category selection (number, name, or prompt sentinel)
    
    Returns:
        The resolved category name
    
    Raises:
        typer.Exit: If invalid category number or no TTY for interactive mode
    """
    from seqnado.tools import get_categories

    categories = get_categories()

    if category == "__PROMPT__":
        if not sys.stdin.isatty():
            logger.error("Interactive category selection requires a TTY.")
            raise typer.Exit(code=1)

        if _RICH_CONSOLE:
            _RICH_CONSOLE.print(
                "[bold cyan]Available Tool Categories:[/bold cyan]"
            )
            for i, cat in enumerate(categories, 1):
                _RICH_CONSOLE.print(f"  {i}. {cat}")
        else:
            typer.echo("Available Tool Categories:")
            for i, cat in enumerate(categories, 1):
                typer.echo(f"  {i}. {cat}")

        category = typer.prompt("Select category number or name")

    if category.isdigit():
        idx = int(category) - 1
        if 0 <= idx < len(categories):
            return categories[idx]
        logger.error(f"Invalid category number: {category} (1-{len(categories)})")
        raise typer.Exit(code=1)
    return category


class _OptionalCategoryCommand(typer.core.TyperCommand):
    """Allow ``--category``/``-c`` to be used with or without a value.

    When the option appears without a following value (end of argv or next
    token starts with ``-``), the sentinel ``__PROMPT__`` is injected so that
    Click sees a valid option-value pair and the command body can trigger an
    interactive prompt.
    """

    def parse_args(self, ctx: click.Context, args: list) -> list:
        """Parse arguments, injecting sentinel for --category without value."""
        patched: list = []
        i = 0
        while i < len(args):
            patched.append(args[i])
            if args[i] in ("--category", "-c"):
                next_val = args[i + 1] if i + 1 < len(args) else None
                if next_val is None or next_val.startswith("-"):
                    patched.append("__PROMPT__")
            i += 1
        return super().parse_args(ctx, patched)
