from pydantic import BaseModel, Field
import pandas as pd
from seqnado.utils import warn_once

class SampleGroup(BaseModel):
    """A single group of samples with an optional reference sample."""
    name: str
    samples: list[str]
    reference_sample: str | None = None

    def __len__(self): return len(self.samples)
    def __contains__(self, sample): return sample in self.samples
    def __str__(self): return f"{self.name}: {len(self.samples)} samples (ref={self.reference_sample})"


class SampleGroups(BaseModel):
    """
    A collection of SampleGroup instances that represent one grouping scheme
    (e.g., for normalization or scaling).
    """
    groups: list[SampleGroup]

    @classmethod
    def from_dataframe(
        cls,
        df: pd.DataFrame,
        subset_column: str = "scaling_group",
        *,
        reference_sample: str | None = None,
    ) -> "SampleGroups":
        """
        Build multiple SampleGroups from a DataFrame based on a grouping column.
        """
        if subset_column not in df.columns:
            warn_once(f"Column '{subset_column}' not found in DataFrame. Returning empty SampleGroups.")
            return cls(groups=[])

        groups = []
        for group_value, group_df in df.groupby(subset_column):
            sample_names = group_df.index.tolist()
            ref_sample = reference_sample if reference_sample in sample_names else sample_names[0]
            groups.append(SampleGroup(name=str(group_value), samples=sample_names, reference_sample=ref_sample))

        return cls(groups=groups)
    
    @property
    def group_names(self) -> list[str]:
        return [g.name for g in self.groups]

    @property
    def sample_to_group(self) -> dict[str, str]:
        return {sample: g.name for g in self.groups for sample in g.samples}

    @property
    def group_to_samples(self) -> dict[str, list[str]]:
        return {g.name: g.samples for g in self.groups}

    def get_group(self, name: str) -> SampleGroup:
        for group in self.groups:
            if group.name == name:
                return group
        raise KeyError(f"Group '{name}' not found.")

    def get_samples(self, group_name: str) -> list[str]:
        return self.get_group(group_name).samples
    
    @property
    def empty(self) -> bool:
        """Check if there are no groups defined."""
        return len(self.groups) == 0

    def __str__(self):
        return f"SampleGroups({len(self.groups)} groups: {', '.join(g.name for g in self.groups)})"
    
    def __len__(self):
        return len(self.groups)


class SampleGroupings(BaseModel):
    """
    A container for multiple named SampleGroups sets,
    e.g., for 'normalization', 'scaling', 'visualization', etc.
    """
    groupings: dict[str, SampleGroups] = Field(default_factory=dict)

    def add_grouping(self, name: str, groups: SampleGroups):
        self.groupings[name] = groups

    def get_grouping(self, name: str) -> SampleGroups:
        if name not in self.groupings:
            raise KeyError(f"Grouping '{name}' not found.")
        return self.groupings[name]

    def __contains__(self, name: str) -> bool:
        return name in self.groupings

    def __str__(self):
        return f"SampleGroupings({', '.join(self.groupings)})"

    def all_samples(self) -> list[str]:
        """Return all unique samples across all groupings."""
        return list({sample for g in self.groupings.values() for group in g.groups for sample in group.samples})
    
    @property
    def empty(self) -> bool:
        """Check if there are no groupings defined."""
        return len(self.groupings) == 0
