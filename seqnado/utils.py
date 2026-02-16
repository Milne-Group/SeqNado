"""Utility functions for SeqNado CLI and general operations."""
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Pattern, Tuple, Union

from loguru import logger

# Compatible typing for compiled regex across Python versions
try:
    RegexType = re.Pattern  # type: ignore[name-defined]
except AttributeError:
    RegexType = type(re.compile(""))

def extract_cores_from_options(options: List[str]) -> Tuple[List[str], int]:
    """
    Extract the number of cores from the snakemake options.
    """

    try:
        cores_flag = options.index("-c")
        cores = int(options[cores_flag + 1])
        options = [
            o for i, o in enumerate(options) if i not in [cores_flag, cores_flag + 1]
        ]
    except ValueError:
        try:
            cores_flag = options.index("--cores")
            cores = int(options[cores_flag + 1])
            options = [
                o
                for i, o in enumerate(options)
                if i not in [cores_flag, cores_flag + 1]
            ]
        except ValueError:
            cores = 1
    except IndexError:
        cores = 1
        options = [o for i, o in enumerate(options) if i not in [cores_flag]]
        logger.warning("Core flag provided but no value given. Defaulting to 1 core.")

    return options, cores


def extract_apptainer_args(options: List[str]) -> Tuple[List[str], str]:
    """
    Extract the apptainer arguments from the snakemake options.
    """
    try:
        apptainer_flag = options.index("--apptainer-args")
        apptainer_args = options[apptainer_flag + 1]
        options = [
            o
            for i, o in enumerate(options)
            if i not in [apptainer_flag, apptainer_flag + 1]
        ]
    except ValueError:
        apptainer_args = ""

    return options, apptainer_args


def remove_unwanted_run_files():
    import glob
    import os
    import shutil

    slurm_files = glob.glob("slurm-*.out")
    sps_files = glob.glob("sps-*")
    simg_files = glob.glob("*.simg")

    for fn in [*slurm_files, *sps_files, *simg_files]:
        try:
            if not os.path.isdir(fn):
                os.remove(fn)
            else:
                shutil.rmtree(fn)

        except Exception as e:
            print(e)


def run_batch_job_on_error(email: str):
    """
    Run a batch job with slurm to send an email on error.
    """

    import subprocess

    slurm_script = f"""#!/bin/bash
    #SBATCH --job-name=seqnado_error_notification
    #SBATCH --mail-type=END
    #SBATCH --mail-user={email}

    echo "An error occurred in the job. Please check the logs for more details."
    """

    with open("error_notification.sh", "w") as f:
        f.write(slurm_script)

    try:
        subprocess.run(["sbatch", "error_notification.sh"], check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to submit slurm job: {e}")


def pepe_silvia():
    print("PEPE SILVIA")
    _pepe_silvia = "https://relix.com/wp-content/uploads/2017/03/tumblr_o16n2kBlpX1ta3qyvo1_1280.jpg"
    return _pepe_silvia



def _as_list(
    items: Union[str, Iterable[Union[str, Pattern]], Pattern, None]
) -> List[Union[str, Pattern]]:
    """Normalize possible inputs into a list of strings or compiled patterns."""
    if items is None:
        return []
    # If it's a compiled pattern or a single string, wrap in list
    if isinstance(items, re.Pattern) or isinstance(items, str):
        return [items]
    # Otherwise assume it's iterable of strings/patterns
    return list(items)


class FileSelector:
    """
    Encapsulates a collection of file names and provides flexible selection.

    Main method:
        select(
            suffix: str,
            includes=None,
            excludes=None,
            case_sensitive=False,
            includes_all=False,
            use_regex=False,
        ) -> List[str]

    - includes / excludes accept strings, iterables of strings, or compiled re.Pattern objects.
    - If use_regex=True string items are treated as regex patterns (compiled with IGNORECASE if case_sensitive=False).
    - If use_regex=False string items are treated as simple substrings (case sensitivity controlled by case_sensitive).
    - includes_all=True requires ALL include descriptors to match (AND). If False, ANY match suffices (OR).
    """

    def __init__(self, files: Iterable[str] = ()):
        self.files: List[str] = list(files)

    @classmethod
    def from_path(cls, path: Union[str, Path], *, recursive: bool = False) -> "FileSelector":
        """
        Create a FileSelector by listing filenames under `path`.
        - If path is a file, returns selector with just that file path (string).
        - If path is a directory, lists child file paths (full paths as strings). If recursive=True, walks recursively.
        """
        p = Path(path)
        if p.is_file():
            return cls([str(p)])
        if not p.exists():
            return cls([])

        if recursive:
            files = [str(fp) for fp in p.rglob("*") if fp.is_file()]
        else:
            files = [str(fp) for fp in p.iterdir() if fp.is_file()]

        return cls(files)

    def add(self, *files: str) -> None:
        """Append files to the selector's pool."""
        self.files.extend(files)

    def clear(self) -> None:
        """Remove all files."""
        self.files.clear()

    # ---------- Matching helpers (top-level methods, no nested functions) ----------
    def _preprocess_items(
        self,
        items: Union[str, Iterable[Union[str, Pattern]], Pattern, None],
        case_sensitive: bool,
        use_regex: bool,
    ) -> List[Dict[str, Any]]:
        """
        Convert user-provided items into a list of processed descriptors:
        - {'kind': 'regex', 'pattern': compiled_pattern}
        - {'kind': 'substr', 'text': processed_substring}
        """
        raw = _as_list(items)
        processed: List[Dict[str, Any]] = []

        for it in raw:
            # If user passed a compiled regex
            if isinstance(it, re.Pattern):
                patt = it
                # Ensure case-insensitive flag if requested
                if not case_sensitive and not (patt.flags & re.IGNORECASE):
                    patt = re.compile(patt.pattern, patt.flags | re.IGNORECASE)
                processed.append({"kind": "regex", "pattern": patt})
                continue

            # plain string
            assert isinstance(it, str)
            if use_regex:
                flags = 0 if case_sensitive else re.IGNORECASE
                patt = re.compile(it, flags)
                processed.append({"kind": "regex", "pattern": patt})
            else:
                # simple substring matching; store normalized substring
                text = it if case_sensitive else it.lower()
                processed.append({"kind": "substr", "text": text})

        return processed

    def _matches_descriptor(self, descriptor: Dict[str, Any], fname: str, fname_for_substr: str) -> bool:
        """
        Test whether a single processed descriptor matches the file name.
        - For 'regex': use descriptor['pattern'].search on original fname.
        - For 'substr': test descriptor['text'] in fname_for_substr.
        """
        if descriptor["kind"] == "regex":
            patt: Pattern = descriptor["pattern"]
            return bool(patt.search(fname))
        else:
            text: str = descriptor["text"]
            return text in fname_for_substr

    # ---------- Public API ----------
    def select(
        self,
        suffix: str,
        includes: Union[str, Iterable[Union[str, Pattern]], Pattern, None] = None,
        excludes: Union[str, Iterable[Union[str, Pattern]], Pattern, None] = None,
        *,
        case_sensitive: bool = False,
        includes_all: bool = False,
        use_regex: bool = False,
    ) -> List[str]:
        """
        Filter files by suffix, include/exclude lists, optional regex, and AND/OR semantics.

        Args:
            suffix: file suffix to match (e.g. ".txt" or "csv"). Leading dot optional.
            includes/excludes: substrings or regex patterns (or compiled patterns).
            case_sensitive: whether substring or regex matching is case-sensitive.
            includes_all: if True, require ALL include descriptors match (AND). If False, any match suffices (OR).
            use_regex: if True, treat string includes/excludes as regex patterns.

        Returns:
            List[str]: filtered list of filenames (the original strings from self.files).
        """
        # normalize suffix: keep leading dot optional, compare on lower() if not case_sensitive
        if not suffix:
            suffix_proc = ""
        else:
            # ensure leading dot if user provided 'txt' vs '.txt' — we'll treat both the same
            if not suffix.startswith("."):
                suffix = "." + suffix
            suffix_proc = suffix if case_sensitive else suffix.lower()

        includes_proc = self._preprocess_items(includes, case_sensitive, use_regex)
        excludes_proc = self._preprocess_items(excludes, case_sensitive, use_regex)

        result: List[str] = []
        for f in self.files:
            fname = f  # original filename for regex matching
            fname_for_substr = f if case_sensitive else f.lower()

            # suffix check
            if suffix_proc:
                if not fname_for_substr.endswith(suffix_proc):
                    continue

            # includes check
            if includes_proc:
                if includes_all:
                    if not all(self._matches_descriptor(d, fname, fname_for_substr) for d in includes_proc):
                        continue
                else:
                    if not any(self._matches_descriptor(d, fname, fname_for_substr) for d in includes_proc):
                        continue

            # excludes check (none of them may match)
            if excludes_proc:
                if any(self._matches_descriptor(d, fname, fname_for_substr) for d in excludes_proc):
                    continue

            result.append(f)

        return result

# ----------------------- CLI Profile & Flag Helpers -----------------------


def create_flag_filter(allowed_flags: Tuple[str, ...]):
    """
    Factory function to create a flag filtering function.
    
    Args:
        allowed_flags: Tuple of allowed flag names (e.g., ("-n", "--dry-run", ...))
    
    Returns:
        A function that checks if an option matches allowed flags.
    """
    def should_pass_flag(opt: str) -> bool:
        """Check if option matches any of the allowed flags."""
        for p in allowed_flags:
            if opt == p or opt.startswith(p + "="):
                return True
        return False
    return should_pass_flag


def get_profile_name(fn: Path) -> Optional[str]:
    """
    Extract profile shortcode from profile directory name.
    
    E.g., "profile_local_conda" -> "lc", "profile_slurm_singularity" -> "ss"
    
    Args:
        fn: Path to profile directory
    
    Returns:
        Profile shortcode or None if not a valid profile name
    """
    name = fn.name
    if not name.startswith("profile_"):
        return None

    profile_parts = name.split("_")
    if len(profile_parts) < 2:
        return None
    initials = "".join(part[0] for part in profile_parts[1:] if part)
    return initials


def get_preset_profiles() -> Dict[str, str]:
    """
    Discover and map all available Snakemake profile presets.
    
    Returns:
        Dict mapping shortcode (e.g., "lc") to profile directory name
    """
    from importlib import resources
    
    profiles_trav = resources.files("seqnado.workflow.envs.profiles")
    profiles = [
        f.name
        for f in profiles_trav.iterdir()
        if f.is_dir()
    ]

    # Map profile shortcuts to directory names
    # E.g., "profile_local_conda" -> "lc", "profile_slurm_singularity" -> "ss"
    return {
        get_profile_name(Path(p)): p 
        for p in profiles 
        if p.startswith("profile_") and get_profile_name(Path(p))
    }


def resolve_profile_path(
    preset: Optional[str],
    pkg_root_trav=None,  # Optional importlib.resources Traversable
):
    """
    Resolve a profile preset to a path.
    
    Prefers profiles installed during `seqnado init` over bundled profiles.
    
    Checks in order:
    1. User's ~/.config/snakemake/{profile_dir_name} (installed via `seqnado init`)
    2. Bundled profile in the package (via pkg_root_trav)
    
    Args:
        preset: Profile preset shortcode (e.g., "lc", "le", "ss", None)
        pkg_root_trav: Optional package traversable (seqnado package)
    
    Returns:
        Path (for user-installed) or Traversable (for bundled), or None if not found
    """
    if not preset:
        return None
    
    profiles = get_preset_profiles()
    profile_dir_name = profiles.get(preset.lower())
    
    if not profile_dir_name:
        return None
    
    # First, check user's ~/.config/snakemake/ for profiles installed via `seqnado init`
    user_profile_path = Path.home() / ".config" / "snakemake" / profile_dir_name
    if user_profile_path.exists():
        return user_profile_path
    
    # Fall back to bundled profiles if pkg_root_trav is provided
    if pkg_root_trav:
        return (
            pkg_root_trav.joinpath("workflow")
            .joinpath("envs")
            .joinpath("profiles")
            .joinpath(profile_dir_name)
        )
    
    return None