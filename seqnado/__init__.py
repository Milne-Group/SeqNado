from .core import (
    Assay,
    NONE_VALUES,
    DataScalingTechnique,
    PileupMethod,
    PCRDuplicateHandling,
    PCRDuplicateTool,
    PeakCallingMethod,
    MotifMethod,
    SpikeInMethod,
    SNPCallingMethod,
    QuantificationMethod,
    MethylationMethod,
    AssaysWithPeakCalling,
    AssaysWithHeatmaps,
    AssaysWithSpikein,
    Molecule,
    Organism,
    FileType,
)

_LAZY_SUBMODULES = {"data", "config", "inputs", "outputs", "analysis"}
_LAZY_ATTRS = {
    "GenomicCoordinate": ("config", "GenomicCoordinate"),
    "SeqNadoProject": ("analysis", "SeqNadoProject"),
    "SeqNadoMultiProject": ("analysis", "SeqNadoMultiProject"),
    "open_project": ("analysis", "open_project"),
}


def __getattr__(name: str):
    if name in _LAZY_SUBMODULES:
        import importlib
        mod = importlib.import_module(f".{name}", __name__)
        globals()[name] = mod
        return mod
    if name in _LAZY_ATTRS:
        mod_name, attr = _LAZY_ATTRS[name]
        import importlib
        mod = importlib.import_module(f".{mod_name}", __name__)
        obj = getattr(mod, attr)
        globals()[name] = obj
        return obj
    raise AttributeError(f"module 'seqnado' has no attribute {name!r}")


__all__ = [
    "Assay",
    "ILLUMINA_FILENAME_PATTERNS",
    "NONE_VALUES",
    "DataScalingTechnique",
    "PileupMethod",
    "INPUT_CONTROL_SUBSTRINGS",
    "PCRDuplicateHandling",
    "PCRDuplicateTool",
    "PeakCallingMethod",
    "MotifMethod",
    "SpikeInMethod",
    "SNPCallingMethod",
    "QuantificationMethod",
    "MethylationMethod",
    "SeqNadoProject",
    "SeqNadoMultiProject",
    "open_project",
    "analysis",
    "data",
    "config",
    "inputs",
    "outputs",
    "AssaysWithPeakCalling",
    "AssaysWithHeatmaps",
    "AssaysWithSpikein",
    "Molecule",
    "Organism",
    "FileType",
    "GenomicCoordinate",
]