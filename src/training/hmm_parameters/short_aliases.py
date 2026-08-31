from src.training.hmm_parameters.training_types import HMM_Chosen_Steps_Schema
from src.antismash_parsing.genomic_context import (
    ModuleGenomicContextFeature,
    ModuleGenomicContext,
)
from src.hmm.hmm_auxiliary_types import (
    DetailedHMMEdgeType,
    DetailedHMMStateType
)

MCF = ModuleGenomicContextFeature
GC = ModuleGenomicContext
ST = DetailedHMMStateType
ET = DetailedHMMEdgeType
SC = HMM_Chosen_Steps_Schema
