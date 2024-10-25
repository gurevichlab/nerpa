@dataclass
class MonomerInsertStepInfo:
    nrp_id: str
    dist_to_prev_module: int  # is negative for insertions in the very beginning of the BGC
    prev_module_context: ModuleLocFeatures  # or the first module in the BGC if dist_to_prev_module is negative
    prev_module_info: AlignmentStep_BGC_Module_Info  # or the first module in the BGC if dist_to_prev_module is negative

    def __eq__(self, other):
        return all([
            self.prev_module_info.gene_id == other.prev_module_info.gene_id,
            self.prev_module_info.a_domain_idx == other.prev_module_info.a_domain_idx,
            self.dist_to_prev_module == other.dist_to_prev_module,
            ])

