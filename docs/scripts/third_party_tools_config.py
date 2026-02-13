"""
Dynamic documentation functions for ThirdPartyToolsConfig.

These functions can be called during MkDocs build process to generate
up-to-date documentation from the current model definitions.
"""

from typing import Dict, List, Any, Type, get_origin, get_args, Optional, Union
from enum import Enum
import inspect
import textwrap

from pydantic import BaseModel
from pydantic.fields import FieldInfo


def get_tool_categories() -> Dict[str, List[str]]:
    """
    Define tool categories for documentation organization.
    Update this function when adding new tools or reorganizing categories.
    """
    return {
        "Alignment Tools": ["bowtie2", "star"],
        "Processing Tools": ["samtools", "picard"], 
        "Trimming Tools": ["cutadapt", "trimgalore"],
        "Peak Calling Tools": ["macs", "lanceotron", "lanceotronmcc", "seacr"],
        "Analysis Tools": ["deeptools", "homer", "bamnado", "subread"],
        "Specialized Tools": ["methyldackel", "bcftools"]
    }


def extract_field_type(field_info: FieldInfo) -> Optional[Type]:
    """Extract the actual type from a field, handling Optional and other generics."""
    annotation = field_info.annotation
    
    if hasattr(annotation, '__origin__'):
        # Handle Union types (like Optional[T] which is Union[T, None])
        if get_origin(annotation) is Union:
            args = get_args(annotation)
            # Find the non-None type
            non_none_args = [arg for arg in args if arg is not type(None)]
            if non_none_args:
                return non_none_args[0]
    
    return annotation


def get_default_value_info(field_info: FieldInfo) -> Dict[str, Any]:
    """Extract default value information from a field."""
    info = {
        "has_default": field_info.default is not None,
        "default_value": field_info.default,
        "has_default_factory": field_info.default_factory is not None,
        "is_required": field_info.is_required()
    }
    
    # Try to get default factory result for documentation
    if field_info.default_factory:
        try:
            info["factory_example"] = field_info.default_factory()
        except Exception:
            info["factory_example"] = None
    
    return info


def analyze_tool_config(tool_class: Type[BaseModel]) -> Dict[str, Any]:
    """Analyze a tool configuration class and extract documentation info."""
    info = {
        "name": tool_class.__name__,
        "docstring": inspect.getdoc(tool_class) or "",
        "fields": {},
        "is_tool_config": issubclass(tool_class, BaseModel),
        "inheritance": [base.__name__ for base in tool_class.__bases__],
        "enums": {},
        "sub_configs": {}
    }
    
    # Analyze fields
    for field_name, field_info in tool_class.model_fields.items():
        field_type = extract_field_type(field_info)
        default_info = get_default_value_info(field_info)
        
        field_analysis = {
            "description": field_info.description or "",
            "type": str(field_type) if field_type else "Unknown",
            "default_info": default_info
        }
        
        # Check if it's an enum
        if field_type and inspect.isclass(field_type) and issubclass(field_type, Enum):
            field_analysis["enum_values"] = [e.value for e in field_type]
            info["enums"][field_type.__name__] = {
                "values": [(e.name, e.value) for e in field_type],
                "docstring": inspect.getdoc(field_type) or ""
            }
        
        # Check if it's a sub-configuration
        if field_type and inspect.isclass(field_type) and issubclass(field_type, BaseModel):
            field_analysis["is_sub_config"] = True
            info["sub_configs"][field_name] = analyze_tool_config(field_type)
        
        info["fields"][field_name] = field_analysis
    
    return info


def generate_tool_documentation_markdown(tool_name: str, tool_class: Type[BaseModel], field_info: FieldInfo) -> str:
    """Generate markdown documentation for a single tool."""
    analysis = analyze_tool_config(tool_class)
    
    md = f"""# {tool_name.title()}

**{field_info.description or 'Tool configuration'}**

{analysis['docstring']}

## Configuration Structure

"""
    
    # Add inheritance info if relevant
    if analysis['inheritance'] and analysis['inheritance'] != ['BaseModel']:
        md += f"**Inherits from:** {', '.join(analysis['inheritance'])}\n\n"
    
    # Document main fields
    if analysis['fields']:
        md += "### Fields\n\n"
        for field_name, field_data in analysis['fields'].items():
            md += f"#### `{field_name}`\n\n"
            
            if field_data.get('description'):
                md += f"{field_data['description']}\n\n"
            
            md += f"**Type:** `{field_data['type']}`\n\n"
            
            # Default value information
            default_info = field_data['default_info']
            if default_info['is_required']:
                md += "**Required:** Yes\n\n"
            elif default_info['has_default']:
                md += f"**Default:** `{default_info['default_value']}`\n\n"
            elif default_info['has_default_factory']:
                if default_info.get('factory_example'):
                    md += f"**Default:** Generated by factory function\n\n"
                else:
                    md += f"**Default:** Factory-generated\n\n"
            
            # Enum values
            if 'enum_values' in field_data:
                md += f"**Possible values:** {', '.join(f'`{v}`' for v in field_data['enum_values'])}\n\n"
            
            # Sub-configuration
            if field_data.get('is_sub_config'):
                md += f"**Sub-configuration:** See [{field_name} details](#{field_name.lower()}-configuration)\n\n"
    
    # Document enums
    if analysis['enums']:
        md += "## Enums\n\n"
        for enum_name, enum_data in analysis['enums'].items():
            md += f"### {enum_name}\n\n"
            if enum_data['docstring']:
                md += f"{enum_data['docstring']}\n\n"
            
            md += "| Name | Value |\n|------|-------|\n"
            for name, value in enum_data['values']:
                md += f"| `{name}` | `{value}` |\n"
            md += "\n"
    
    # Document sub-configurations
    if analysis['sub_configs']:
        md += "## Sub-configurations\n\n"
        for sub_name, sub_analysis in analysis['sub_configs'].items():
            md += f"### {sub_name.title()} Configuration\n\n"
            
            if sub_analysis['docstring']:
                md += f"{sub_analysis['docstring']}\n\n"
            
            if sub_analysis['fields']:
                md += "#### Fields\n\n"
                for field_name, field_data in sub_analysis['fields'].items():
                    md += f"**`{field_name}`**"
                    if field_data.get('description'):
                        md += f": {field_data['description']}"
                    
                    default_info = field_data['default_info']
                    if not default_info['is_required']:
                        if default_info['has_default']:
                            md += f" (default: `{default_info['default_value']}`)"
                        elif default_info['has_default_factory']:
                            md += f" (default: factory-generated)"
                    
                    md += "\n\n"
    
    return md


def generate_tools_index_markdown(config_class: Type[BaseModel]) -> str:
    """Generate the tools index page with categorized tool list."""
    categories = get_tool_categories()
    
    md = """# Tool Reference

This page provides comprehensive documentation for all bioinformatics tools supported by the configuration system.

## Overview

Each tool configuration follows consistent patterns:
- **threads**: Number of CPU threads to allocate  
- **options**: CLI options string with validation and filtering
- **sub-configurations**: Some tools have multiple operational modes

## Tool Categories

"""
    
    # Get all configured tools
    all_tools = set()
    for field_name, field_info in config_class.model_fields.items():
        tool_class = extract_field_type(field_info)
        if tool_class and inspect.isclass(tool_class) and issubclass(tool_class, BaseModel):
            all_tools.add(field_name)
    
    # Generate categorized listings
    for category, tools in categories.items():
        md += f"### {category}\n\n"
        
        for tool in tools:
            if tool in all_tools:
                field_info = config_class.model_fields[tool]
                description = field_info.description or "Tool configuration"
                md += f"- **[{tool}](tools.md#{tool})** - {description}\n"
        
        md += "\n"
    
    # List any uncategorized tools
    categorized_tools = set()
    for tools in categories.values():
        categorized_tools.update(tools)
    
    uncategorized = all_tools - categorized_tools
    if uncategorized:
        md += "### Other Tools\n\n"
        for tool in sorted(uncategorized):
            field_info = config_class.model_fields[tool]
            description = field_info.description or "Tool configuration"
            md += f"- **[{tool}](tools.md#{tool})** - {description}\n"
        md += "\n"
    
    return md


def generate_complete_tools_documentation(config_class: Type[BaseModel]) -> str:
    """Generate complete documentation for all tools in a single markdown file."""
    md = generate_tools_index_markdown(config_class)
    
    md += """
---

## Detailed Tool Documentation

"""
    
    # Generate documentation for each tool
    for field_name, field_info in config_class.model_fields.items():
        tool_class = extract_field_type(field_info)
        if tool_class and inspect.isclass(tool_class) and issubclass(tool_class, BaseModel):
            md += f"\n---\n\n"
            tool_md = generate_tool_documentation_markdown(field_name, tool_class, field_info)
            # Add anchor for the tool
            tool_md = tool_md.replace(f"# {field_name.title()}", f"## {field_name} {{#{field_name}}}")
            md += tool_md
    
    return md


def generate_assay_documentation(get_assay_specific_tools_func, assay_enum_class) -> str:
    """Generate documentation showing which tools are used for each assay type."""
    md = """# Assay-Specific Tool Configurations

This page shows which tools are automatically configured for each assay type when using the `for_assay()` factory method.

## Assay Types

"""
    
    for assay in assay_enum_class:
        tools = get_assay_specific_tools_func(assay)
        tool_names = [tool_class.__name__.lower() for tool_class in tools]
        
        md += f"### {assay.value.upper()}\n\n"
        md += f"**Tools configured:** {', '.join(f'`{name}`' for name in sorted(tool_names))}\n\n"
        
        # Add brief description of what this assay type is for
        assay_descriptions = {
            "atac": "Assay for Transposase-Accessible Chromatin - identifies open chromatin regions",
            "chip": "Chromatin Immunoprecipitation sequencing - maps protein-DNA interactions", 
            "cat": "Cleavage Under Targets and Tagmentation - low-input chromatin profiling",
            "rna": "RNA sequencing - measures gene expression levels",
            "snp": "Single Nucleotide Polymorphism analysis - identifies genetic variants",
            "meth": "DNA methylation analysis - maps methylation patterns",
            "crispr": "CRISPR screening analysis - analyzes guide RNA performance"
        }
        
        if assay.value.lower() in assay_descriptions:
            md += f"{assay_descriptions[assay.value.lower()]}\n\n"
    
    return md


def generate_usage_examples(config_class: Type[BaseModel]) -> str:
    """Generate usage examples documentation."""
    class_name = config_class.__name__
    
    md = f"""# Usage Examples

## Basic Usage

### Create Configuration for Specific Assay

```python
from seqnado import Assay
from your_module import {class_name}

# Create configuration with assay-specific defaults
config = {class_name}.for_assay(Assay.ATAC)

# Export to dictionary
config_dict = config.model_dump(exclude_none=True)
```

### Custom Tool Configuration

```python
# Override specific tools during creation
config = {class_name}.for_assay(
    Assay.ATAC,
    bowtie2=Bowtie2(
        align=ToolConfig(threads=16, options="--very-fast")
    ),
    samtools=Samtools(
        sort=ToolConfig(threads=12, options="-@ {{threads}} -m 4G")
    )
)
```

### Manual Configuration

```python
# Create empty configuration and add tools manually
config = {class_name}()

# Add specific tools
config.bowtie2 = Bowtie2()
config.samtools = Samtools()
config.deeptools = Deeptools()

# Check which tools are configured
configured_tools = config.get_configured_tools()
print(f"Configured tools: {{list(configured_tools.keys())}}")
```

## Advanced Usage

### Option Filtering

```python
from your_module import OptionsBase, ToolConfig

# Create options with exclusions
options = OptionsBase(
    value="--very-sensitive --threads 8 --fast",
    exclude={{"--fast"}}  # This option will be filtered out
)

# The filtered result will be: "--very-sensitive --threads 8"
print(options.option_string_filtered)
```

### Dynamic Thread Configuration

```python
# Use template strings for dynamic values
config = {class_name}.for_assay(
    Assay.RNA,
    samtools=Samtools(
        sort=ToolConfig(
            threads=16,
            options="-@ {{threads}} -m 2G"  # {{threads}} will be replaced
        )
    )
)
```

### Serialization Options

```python
# Standard dictionary export
config_dict = config.model_dump(exclude_none=True)

# Include descriptions as help fields
config_with_help = config.dump_with_descriptions()

# Export to YAML with comments (requires PyYAML)
yaml_output = config.dump_yaml_with_comments()

# Export to TOML with comments (requires tomli-w)  
toml_output = config.dump_toml_with_comments()
```

## Integration Examples

### With Configuration Files

```python
import json
from pathlib import Path

# Save configuration
config = {class_name}.for_assay(Assay.CHIP)
config_path = Path("config.json")
config_path.write_text(config.model_dump_json(indent=2, exclude_none=True))

# Load configuration
config_data = json.loads(config_path.read_text())
loaded_config = {class_name}(**config_data)
```

### With Environment-Specific Overrides

```python
import os

def create_config_for_environment():
    # Base configuration
    config = {class_name}.for_assay(Assay.ATAC)
    
    # Override based on environment
    if os.getenv("HIGH_MEMORY") == "true":
        config.samtools.sort.options = "-@ {{threads}} -m 8G"
    
    if os.getenv("FAST_MODE") == "true":
        config.bowtie2.align.align.options = "--very-fast"
    
    return config
```

### Validation and Error Handling

```python
try:
    # This will raise validation error for invalid options
    config = {class_name}(
        bowtie2=Bowtie2(
            align=ToolConfig(
                threads=0,  # Invalid: less than minimum
                options="--invalid'quote"  # Invalid: bad quoting
            )
        )
    )
except ValueError as e:
    print(f"Configuration error: {{e}}")
```
"""
    return md


# Convenience function to generate all documentation
def generate_all_documentation(config_class: Type[BaseModel], 
                             get_assay_specific_tools_func,
                             assay_enum_class) -> Dict[str, str]:
    """
    Generate all documentation pages.
    
    Returns a dictionary mapping filename to content for each documentation page.
    Note: tools-index.md is not generated as docs/tools.md serves as the main reference.
    """
    return {
        "assays.md": generate_assay_documentation(get_assay_specific_tools_func, assay_enum_class),
        "tool-config-examples.md": generate_usage_examples(config_class)
    }


# Example usage for MkDocs integration
def mkdocs_integration_example():
    """
    Example of how to integrate with MkDocs using mkdocs-gen-files plugin.
    
    Add to your mkdocs.yml:
    
    plugins:
      - gen-files:
          scripts:
            - docs/generate_docs.py
    
    Then create docs/generate_docs.py:
    
    ```python
    import mkdocs_gen_files
    from your_module import ThirdPartyToolsConfig, get_assay_specific_tools, Assay
    from docs_generator import generate_all_documentation
    
    # Generate all documentation
    docs = generate_all_documentation(
        ThirdPartyToolsConfig, 
        get_assay_specific_tools, 
        Assay
    )
    
    # Write files
    for filename, content in docs.items():
        with mkdocs_gen_files.open(filename, "w") as f:
            f.write(content)
    ```
    """
    pass
# Generate documentation if running under mkdocs-gen-files
try:
    import mkdocs_gen_files
    from seqnado.config.third_party_tools import ThirdPartyToolsConfig, get_assay_specific_tools
    from seqnado import Assay

    # Generate all documentation
    docs = generate_all_documentation(
        ThirdPartyToolsConfig, 
        get_assay_specific_tools, 
        Assay
    )
    
    # Write files
    for filename, content in docs.items():
        with mkdocs_gen_files.open(filename, "w") as f:
            f.write(content)
except ImportError:
    pass
