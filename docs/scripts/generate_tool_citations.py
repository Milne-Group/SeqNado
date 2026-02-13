#!/usr/bin/env python3
"""
Generate Tool Sections for Documentation

This script uses SeqNado's tools module to automatically generate
tool sections in docs/tools.md, organized by category.

Usage:
    python docs/scripts/generate_tool_citations.py --update
    python docs/scripts/generate_tool_citations.py --dry-run
"""

import argparse
import re
import sys
from pathlib import Path
from typing import Dict

from seqnado.tools import (
    format_citation,
    get_available_tools,
    get_categories,
    get_tool_citation,
    get_tool_version,
    get_tool_version_from_container,
    is_apptainer_available,
)

TOOL_TEMPLATE = """#### {display_name}
**Purpose**: {description}  
**Version**: {version}  
**Usage**: {usage}  
**Reference**: {citation}

"""


def _get_tool_version_for_docs(tool_name: str) -> str:
    """Get tool version for documentation.

    Uses ``conda list --json`` from the container (one call per unique
    container, cached) as the primary source.  Falls back to local
    detection when apptainer is not available.
    """
    _NOT_FOUND = "Version information not available"

    # 1. Try conda metadata from container (bulk, cached per container)
    if is_apptainer_available():
        try:
            version = get_tool_version_from_container(tool_name)
            if version:
                return version
        except Exception as e:
            print(
                f"Warning: Container metadata check failed for '{tool_name}': {e}",
                file=sys.stderr,
            )

    # 2. Fall back to local detection
    try:
        version = get_tool_version(tool_name, use_container=False)
        if version and version != _NOT_FOUND:
            return version
    except Exception as e:
        print(
            f"Warning: Local version check failed for '{tool_name}': {e}",
            file=sys.stderr,
        )

    return "Latest via container"


def generate_tool_entry(tool_name: str, tool_info: Dict) -> str:
    """Generate a markdown entry for a single tool"""
    # Get display name from tool info, fallback to formatted tool name
    display_name = tool_info.get(
        "display_name", tool_name.replace("_", " ").replace("-", " ").title()
    )

    # Get description
    description = tool_info.get("description", "")

    # Get usage description from tool info, fallback to description
    usage = tool_info.get("usage", description)

    # Get version via live detection (container or local)
    version = _get_tool_version_for_docs(tool_name)

    # Get citation
    citation_text = "See documentation"
    citation_key = tool_info.get("citation")
    if citation_key:
        try:
            bibtex = get_tool_citation(tool_name)
            if bibtex:
                formatted = format_citation(bibtex)
                if formatted:
                    citation_text = formatted
        except Exception:
            print(
                f"Warning: Could not get citation for tool '{tool_name}'",
                file=sys.stderr,
            )
            pass

    return TOOL_TEMPLATE.format(
        display_name=display_name,
        description=description,
        version=version,
        usage=usage,
        citation=citation_text,
    )


def generate_tools_section() -> str:
    """Generate the complete Tools section organized by category"""
    tools = get_available_tools()
    categories = get_categories()

    lines = [
        "## Tools",
        "",
        "Tools are organized by category, matching the structure of the `seqnado tools` CLI command. For more information about any tool, use `seqnado tools <toolname>`.",
        "",
        "<!-- AUTO-GENERATED TOOL SECTIONS - DO NOT EDIT MANUALLY -->",
        "<!-- This section is automatically updated by docs/scripts/generate_tool_citations.py -->",
        "<!-- To update, run: python docs/scripts/generate_tool_citations.py --update -->",
        "",
    ]

    # Organize tools by category
    for category in categories:
        category_tools = {
            name: info
            for name, info in tools.items()
            if info.get("category") == category
        }

        if not category_tools:
            continue

        lines.append(f"### {category}")
        lines.append("")

        # Sort tools alphabetically within category
        for tool_name in sorted(category_tools.keys()):
            tool_info = category_tools[tool_name]
            entry = generate_tool_entry(tool_name, tool_info)
            lines.append(entry)

    return "\n".join(lines)


def update_tools_md(doc_path: Path, dry_run: bool = False) -> bool:
    """Update the tools.md file with auto-generated tool sections"""
    if not doc_path.exists():
        print(f"Error: {doc_path} does not exist")
        return False

    with open(doc_path, "r", encoding="utf-8") as f:
        content = f.read()

    # Generate new tools section
    print("Generating tool sections from SeqNado API...")
    tools_section = generate_tools_section()

    # Replace from "## Tools" to end of file
    pattern = r"\n## Tools\n.*$"

    if re.search(pattern, content, re.DOTALL):
        new_content = re.sub(pattern, "\n" + tools_section, content, flags=re.DOTALL)
        print("Replaced existing Tools section")
    else:
        # Append to end
        new_content = content.rstrip() + "\n\n" + tools_section + "\n"
        print("Appended new Tools section")

    if dry_run:
        print("\nDry run - changes not saved. Preview of updates:")
        print("=" * 80)
        lines = new_content.split("\n")
        print("\n".join(lines[-100:]))
        print("=" * 80)
        return True

    with open(doc_path, "w", encoding="utf-8") as f:
        f.write(new_content)

    print(f"Successfully updated {doc_path}")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Generate tool sections in documentation from SeqNado's tools API"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("docs/tools.md"),
        help="Path to the markdown file (default: docs/tools.md)",
    )
    parser.add_argument(
        "--update", action="store_true", help="Update the file in place"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without making changes",
    )

    args = parser.parse_args()

    if not args.update and not args.dry_run:
        print("Please specify --update or --dry-run")
        return 1

    update_tools_md(args.output, dry_run=args.dry_run)
    return 0


if __name__ == "__main__":
    exit(main())
