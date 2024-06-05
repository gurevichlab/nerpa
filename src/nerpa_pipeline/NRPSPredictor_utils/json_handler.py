from typing import Union, Literal, TypedDict, List, Dict, Tuple
from src.nerpa_pipeline.NRPSPredictor_utils.antismash_parser_types import (
BGC_Cluster,
Coords,
Gene,
Module,
A_Domain,
C_DOMAIN,
ModifyingDomain,
ConnectingDomain,
SVM_LEVEL,
SVM_Prediction,
STRAND
)
from parse import parse
from collections import defaultdict

GeneId = str

class A_Domain_Id:
    gene_id: GeneId
    idx: int

    def __init__(self, a_domain_id: str):
        parsed_id = parse('nrpspksdomains_{gene_id}_AMP-binding.{idx:d})', a_domain_id)
        self.gene_id = parsed_id['gene_id']
        self.idx = parsed_id['idx']


# tested on antismash v7
def extract_a_domain_specificity_info(a_domain_data: dict) -> A_Domain:
    svm_dict = a_domain_data['nrpys']

    svm = {}
    for svm_level, svm_level_str in [(SVM_LEVEL.PHYSIOCHEMICAL_CLASS, 'physiochemical_class'),
                                     (SVM_LEVEL.LARGE_CLUSTER, 'large_cluster'),
                                     (SVM_LEVEL.SMALL_CLUSTER, 'small_cluster'),
                                     (SVM_LEVEL.SINGLE_AMINO, 'single_amino')]:
        score = svm_dict[svm_level_str]['score']
        substrates = svm_dict[svm_level_str]['substrates']
        svm[svm_level] = SVM_Prediction(score=svm_dict[svm_level_str]['score'],
                                        monomer_residues=[substrate['short'].lower() for substrate in substrates])

    return A_Domain(aa10=svm_dict['aa10'],
                    aa34=svm_dict['aa34'],
                    svm=svm)


def extract_a_domains_info(contig_data: dict) -> Dict[GeneId, List[A_Domain]]:
    if "antismash.modules.nrps_pks" not in contig_data["modules"]:
        return defaultdict(list)

    def extract_a_domain_info(domain_id: str, antismash_prediction: dict) -> Tuple[A_Domain_Id, A_Domain]:
        try:
            parsed_id = A_Domain_Id(domain_id)
            a_domain = extract_a_domain_specificity_info(antismash_prediction)
        except KeyError:
            raise RuntimeError('Unable to parse A-domain prediction. Probably, an old version of antismash is used.')
        return parsed_id, a_domain

    a_domains_with_ids =  [extract_a_domain_info(domain_id, prediction)
                            for domain_id, prediction in contig_data['modules']['antismash.modules.nrps_pks']['domain_predictions'].items()
                            if 'AMP-binding' in domain_id]
    a_domains_per_gene = defaultdict(list)
    for a_domain_id, a_domain in sorted(a_domains_with_ids,
                                        key=lambda a_domain_with_id: a_domain_with_id[0].idx):
        a_domains_per_gene[a_domain_id.gene_id].append(a_domain)
    return a_domains_per_gene


def parse_cds_coordinates(location: str) -> Coords:
    def parse_location(loc: str) -> Coords:
        # e.g. 'loc' = '[351:486](+)'
        parsed_loc = parse('[{start:d}:{end:d}]({strand})', loc)
        return Coords(start=parsed_loc['start'], end=parsed_loc['end'],
                      strand=STRAND.FORWARD if parsed_loc['strand'] == '+' else STRAND.REVERSE)

    if not location.startswith('join{'):
        return parse_location(location)

    # split gene (fungal insertion or whatever)
    location_parts = location[len('join{'):-1].split(',')
    parsed_location_parts = sorted(map(parse_location, location_parts),
                                   key=lambda x: x.start)

    return Coords(start=parsed_location_parts[0].start,  # return merged coordinates
                  end=parsed_location_parts[-1].end,
                  strand=parsed_location_parts[0].strand)  # it is assumed that all parts have the same strand


def extract_gene_coords(contig_data: dict) -> Dict[GeneId, Coords]:
    gene_coords = {}
    for feature in contig_data['features']:
        if feature['type'] == 'gene':
            gene_id = feature['qualifiers']['locus_tag'][0]
            coords = parse_cds_coordinates(feature['location'])
            gene_coords[gene_id] = coords
    return gene_coords


def extract_modules(gene_data: dict, a_domains: List[A_Domain]) -> List[Module]:
    modules = []
    a_domains_iter = iter(a_domains)
    for module_data in gene_data['modules']:
        module = Module()
        for domain_data in module_data['components']:
            match domain_data['hit_id']:
                case 'AMP-binding': module.a_domain = next(a_domains_iter)
                case 'Condensation': module.c_domain = C_DOMAIN.C
                case 'Condensation_Starter': module.c_domain = C_DOMAIN.C_STARTER
                case 'Condensation_LCL': module.c_domain = C_DOMAIN.C_LCL
                case 'Condensation_DCL': module.c_domain = C_DOMAIN.C_DCL
                case 'Condensation_Dual': module.c_domain = C_DOMAIN.C_DUAL
                case 'cMT' | 'nMT' | 'oMT': module.modifying_domains.append(ModifyingDomain.MT)
                case 'Epimerization': module.modifying_domains.append(ModifyingDomain.E)
                case 'Thioesterase' | 'TD': module.terminal_domain = True
                case 'NRPS-COM_Cterm': module.connecting_domain = ConnectingDomain.CTERM
                case 'NRPS-COM_Nterm': module.connecting_domain = ConnectingDomain.NTERM
        modules.append(module)
    return modules


def extract_genes(contig_data: dict, a_domains_per_gene: Dict[GeneId, List[A_Domain]]) -> List[Gene]:
    gene_coords = extract_gene_coords(contig_data)

    genes = []
    for gene_id, gene_data in contig_data['modules']['antismash.detection.nrps_pks_domains']['cds_results'].items():
        modules = extract_modules(gene_data, a_domains_per_gene[gene_id])
        genes.append(Gene(gene_id=gene_id,
                          coords=gene_coords[gene_id],
                          modules=modules))

    genes.sort(key=lambda gene: gene.coords.start)
    return genes


def extract_bgc_clusters(ctg_idx: int, contig_data: dict, genes: List[Gene]) -> List[BGC_Cluster]:
    bgcs = []
    for bgc_data in contig_data['areas']:
        if 'NRPS' in bgc_data['products']:
            bgc_genes = [gene for gene in genes
                         if gene.coords.start <= bgc_data['start'] <= gene.coords.end]
            bgcs.append(BGC_Cluster(contig_id=f'ctg{ctg_idx + 1}',
                                    genes=bgc_genes))
    return bgcs


def parse_antismash_json(antismash_json: dict) -> List[BGC_Cluster]:
    bgcs = []
    for ctg_idx, contig_data in enumerate(antismash_json['records']):
        a_domains_per_gene = extract_a_domains_info(contig_data)
        if not any(a_domains for a_domains in a_domains_per_gene.values()):  # no A-domains found in the contig
            continue

        genes = extract_genes(contig_data, a_domains_per_gene)
        bgcs.extend(extract_bgc_clusters(ctg_idx, contig_data, genes))

    return bgcs