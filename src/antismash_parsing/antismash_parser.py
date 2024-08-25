from typing import List, Dict, Tuple
from src.antismash_parsing.antismash_parser_types import (
    A_Domain,
    antiSMASH_record,
    BGC_Cluster,
    Coords,
    DomainType,
    Gene,
    Module,
    SVM_LEVEL,
    SVM_Prediction,
    STRAND
)
from src.config import antiSMASH_Parsing_Config
from src.monomer_names_helper import antiSMASH_MonomerName
from parse import parse
from collections import defaultdict

GeneId = str


class A_Domain_Id:
    gene_id: GeneId
    idx: int

    def __init__(self, a_domain_id: str):
        parsed_id = parse('nrpspksdomains_{gene_id}_AMP-binding.{idx:d}', a_domain_id)
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
        substrates = svm_dict[svm_level_str]['substrates']
        svm[svm_level] = SVM_Prediction(score=svm_dict[svm_level_str]['score'],
                                        substrates=[antiSMASH_MonomerName(substrate['short'])
                                                    for substrate in substrates])

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

    a_domains_with_ids = [extract_a_domain_info(domain_id, prediction)
                          for domain_id, prediction in contig_data['modules']['antismash.modules.nrps_pks']['domain_predictions'].items()
                          if 'AMP-binding' in domain_id]
    a_domains_per_gene = defaultdict(list)
    for a_domain_id, a_domain in sorted(a_domains_with_ids,
                                        key=lambda a_domain_with_id: a_domain_with_id[0].idx):
        a_domains_per_gene[a_domain_id.gene_id].append(a_domain)
    return a_domains_per_gene


def parse_cds_coordinates(location: str) -> Coords:
    # BGC0000296, gene acmY had location [45371:>47435](+).
    # I have no idea what this is, so I just remove the '<' and '>' and hope for the best
    location = location.replace('<', '').replace('>', '')
    location = location.replace(' ','')  # remove possible spaces

    def parse_location(loc: str) -> Coords:
        # e.g. 'loc' = '[351:486](+)'
        parsed_loc = parse('[{start:d}:{end:d}]({strand})', loc)
        return Coords(start=parsed_loc['start'],
                      end=parsed_loc['end'],
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


def extract_gene_id(feature_qualifiers: dict) -> str:
    for gene_id_key in ['locus_tag', 'gene', 'protein_id']:
        if gene_id_key in feature_qualifiers:
            return feature_qualifiers[gene_id_key][0]
    raise KeyError('Gene ID not found in feature qualifiers')


def extract_gene_coords(contig_data: dict) -> Dict[GeneId, Coords]:
    gene_coords = {}
    for feature in contig_data['features']:
        if feature['type'] in ('gene', 'CDS'):
            try:
                gene_id = extract_gene_id(feature['qualifiers'])
                coords = parse_cds_coordinates(feature['location'])
                gene_coords[gene_id] = coords
            except:
                raise
    return gene_coords


def extract_modules(gene_data: dict, a_domains: List[A_Domain],
                    config: antiSMASH_Parsing_Config) -> List[Module]:
    modules = []
    a_domains_iter = iter(a_domains)
    for module_data in gene_data['modules']:
        module = Module()
        for domain_data in module_data['components']:
            domain_type_str = domain_data['domain']['hit_id']
            if domain_type_str not in config.ANTISMASH_DOMAINS_NAMES:
                continue
            domain_type = DomainType[config.ANTISMASH_DOMAINS_NAMES[domain_type_str]]
            module.domains_sequence.append(domain_type)

            if domain_type == DomainType.A:
                a_domain = next(a_domains_iter)
                module.a_domain = a_domain

        modules.append(module)
    return modules


def extract_genes(contig_data: dict,
                  a_domains_per_gene: Dict[GeneId, List[A_Domain]],
                  config: antiSMASH_Parsing_Config) -> List[Gene]:
    gene_coords = extract_gene_coords(contig_data)

    genes = []
    for gene_id, gene_data in contig_data['modules']['antismash.detection.nrps_pks_domains']['cds_results'].items():
        modules = extract_modules(gene_data, a_domains_per_gene[gene_id], config)
        if modules:
            genes.append(Gene(gene_id=gene_id,
                              coords=gene_coords[gene_id],
                              modules=modules))

    genes.sort(key=lambda gene: gene.coords.start)
    return genes


def extract_bgc_clusters(genome_id: str, ctg_idx: int,
                         contig_data: dict, genes: List[Gene]) -> List[BGC_Cluster]:
    bgcs = []
    for bgc_idx, bgc_data in enumerate(contig_data['areas']):
        if 'NRPS' in bgc_data['products']:
            bgc_genes = [gene for gene in genes
                         if bgc_data['start'] <= gene.coords.start <= gene.coords.end <= bgc_data['end']]

            bgcs.append(BGC_Cluster(genome_id=genome_id,
                                    contig_id=f'ctg{ctg_idx + 1}',
                                    bgc_idx=bgc_idx,
                                    genes=bgc_genes))
    return bgcs


def parse_antismash_json(antismash_json: antiSMASH_record,
                         config: antiSMASH_Parsing_Config) -> List[BGC_Cluster]:
    bgcs = []
    genome_id = antismash_json['input_file']
    for ctg_idx, contig_data in enumerate(antismash_json['records']):
        a_domains_per_gene = extract_a_domains_info(contig_data)
        if not any(a_domains for a_domains in a_domains_per_gene.values()):  # no A-domains found in the contig
            continue

        genes = extract_genes(contig_data, a_domains_per_gene, config)
        bgcs.extend(extract_bgc_clusters(genome_id, ctg_idx, contig_data, genes))

    return bgcs
