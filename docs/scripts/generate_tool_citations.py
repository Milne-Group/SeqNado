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
import subprocess
import sys
from pathlib import Path
from typing import Dict

from seqnado.tools import (
    format_citation,
    get_available_tools,
    get_categories,
    get_tool_citation,
)

TOOL_TEMPLATE = """#### {display_name}
**Purpose**: {description}  
**Version**: {version}  
**Usage**: {usage}  
**Reference**: {citation}

"""


def get_tool_version_from_cli(tool_name: str) -> str:
    """Get tool version by running 'seqnado tools <tool> --version'"""
    try:
        result = subprocess.run(
            ["seqnado", "tools", tool_name, "--version"],
            capture_output=True,
            text=True,
            timeout=60,
        )
        output = result.stdout + result.stderr

        # Remove ANSI color codes
        ansi_escape = re.compile(r"\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])")
        output = ansi_escape.sub("", output)

        # Look for "Version Information:" section
        lines = output.split("\n")
        for i, line in enumerate(lines):
            if "Version Information:" in line:
                if i + 1 < len(lines):
                    version_line = lines[i + 1].strip()

                    # Skip if line is empty or contains help/usage keywords
                    if not version_line or any(
                        skip in version_line.lower()
                        for skip in [
                            "usage:",
                            "command",
                            "help",
                            "option",
                            "--",
                            "positional",
                            "arguments",
                            "subcommand",
                            "for more information",
                        ]
                    ):
                        return "Latest via container"

                    # Return the full version line (contains the version number and any prefix)
                    version_match = re.search(
                        r"(?:v)?\d+(?:\.\d+)+(?:[a-z0-9\-]*)?",
                        version_line,
                        re.IGNORECASE,
                    )
                    if version_match:
                        # Return the full line which contains the version
                        return version_line

                    return "Latest via container"

        return "Latest via container"
    except subprocess.TimeoutExpired:
        return "Latest via container"
    except Exception as e:
        print(
            f"Warning: Could not get version for tool '{tool_name}': {e}",
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

    # Get version by calling seqnado tools command
    version = get_tool_version_from_cli(tool_name)

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
        "Tools are organized by category, matching the structure of the `seqnado tools` CLI command. For more information about any tool, use `seqnado tools <toolname> --info`.",
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
