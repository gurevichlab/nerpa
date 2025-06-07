from __future__ import annotations
from typing import (
    List,
    Dict,
    NamedTuple,
    Tuple,
    Optional,
    Set
)
from src.antismash_parsing.antismash_parser_types import (
    A_Domain,
    antiSMASH_record,
    BGC_Cluster,
    Coords,
    DomainType,
    Gene,
    GeneId,
    Module,
    SVM_LEVEL,
    SVM_Prediction,
    STRAND
)
from src.config import antiSMASH_Processing_Config
from src.data_types import BGC_ID
from src.monomer_names_helper import (
    antiSMASH_MonomerName,
    AA10,
    AA34
)
from parse import parse
from collections import defaultdict
from src.pipeline.logger import NerpaLogger


class A_Domain_Id(NamedTuple):
    gene_id: GeneId
    idx: int

    @classmethod
    def from_antismash_id(cls, a_domain_id: str) -> A_Domain_Id:
        parsed_id = parse('nrpspksdomains_{gene_id}_AMP-binding.{idx:d}', a_domain_id)
        return cls(gene_id=parsed_id['gene_id'], idx=parsed_id['idx'])


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

    return A_Domain(aa10=AA10(svm_dict['aa10']),
                    aa34=AA34(svm_dict['aa34']),
                    svm=svm)


def extract_a_domains_info(contig_data: dict,
                           config: antiSMASH_Processing_Config,
                           ctg_idx: int,  # only for error messages
                           genome_id: str, # only for error messages
                           log: NerpaLogger) -> Dict[A_Domain_Id, A_Domain]:
    if "antismash.modules.nrps_pks" not in contig_data["modules"]:
        return {}

    domain_predictions = contig_data['modules']['antismash.modules.nrps_pks']['domain_predictions']
    a_domain_per_id = {}
    for as_domain_id, antismash_prediction in domain_predictions.items():
        if 'AMP-binding' in as_domain_id:
            try:
                parsed_id = A_Domain_Id.from_antismash_id(as_domain_id)
                specificity_info = extract_a_domain_specificity_info(antismash_prediction)
            except Exception as e:
                if config.DEBUG_MODE:
                    log.error(f'Error parsing A-domain {as_domain_id} '
                              f'in contig {ctg_idx} of genome {genome_id}:\n')
                else:
                    log.warning(f'Error parsing A-domain {as_domain_id} '
                                f'in contig {ctg_idx} of genome {genome_id}')
            else:
                a_domain_per_id[parsed_id] = specificity_info

    return a_domain_per_id


def extract_a_domains_coords(contig_data: dict) -> Dict[A_Domain_Id, Coords]:
    a_domains_coords: Dict[A_Domain_Id, Coords] = {}
    for gene_id, gene_hmm_hits in contig_data['modules']['antismash.detection.nrps_pks_domains']['cds_results'].items():
        gene_a_domains_coords = sorted([Coords.from_hmm_hit(query_start=domain_hmm['query_start'],
                                                            query_end=domain_hmm['query_end'])
                                        for domain_hmm in gene_hmm_hits['domain_hmms']
                                        if domain_hmm['hit_id'] == 'AMP-binding'],
                                       key=lambda x: x.start)
        a_domains_coords.update({A_Domain_Id(gene_id=gene_id, idx=idx): coords
                                 for idx, coords in enumerate(gene_a_domains_coords, start=1)})  # antiSMASH A-domain indices start from 1

    return a_domains_coords


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

    if not any(location.startswith(pref)
                for pref in ['order', 'join']):
        return parse_location(location)

    # split gene (fungal insertion or whatever)
    pref = location[:location.find('{')]  # "order" or "join"
    location_parts = location[len(pref) + 1 : -1].split(',')  # +1 to skip '{' and -1 to skip '}'
    parsed_location_parts = sorted(map(parse_location, location_parts),
                                   key=lambda x: x.start)

    return Coords(start=parsed_location_parts[0].start,  # return merged coordinates
                  end=parsed_location_parts[-1].end,
                  strand=parsed_location_parts[0].strand)  # it is assumed that all parts have the same strand


def extract_gene_id(feature_qualifiers: dict) -> GeneId:
    for gene_id_key in ['locus_tag', 'gene', 'protein_id']:
        if gene_id_key in feature_qualifiers:
            return GeneId(feature_qualifiers[gene_id_key][0])
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


def extract_modules(gene_id: str,  # only for error messages
                    gene_data: dict,
                    a_domains_with_coords: List[Tuple[A_Domain, Coords]],
                    config: antiSMASH_Processing_Config) -> List[Module]:
    def get_domain_sequence(module: dict) -> List[DomainType]:
        return [DomainType[config.ANTISMASH_DOMAINS_NAMES_MAPPING[domain['domain']['hit_id']]]
                for domain in module['components']
                if (domain_name := domain['domain']['hit_id']) in config.ANTISMASH_DOMAINS_NAMES_MAPPING]

    modules = []
    modules_iter = iter(gene_data['modules'])
    a_domains_iter = iter(a_domains_with_coords)
    module, a_domain_with_coords = (next(modules_iter, None), next(a_domains_iter, None))
    while (module is not None) and (a_domain_with_coords is not None):
        a_domain, a_domain_coords = a_domain_with_coords

        if a_domain_coords.end < module['components'][0]['domain']['query_start']:  # orphan A-domain (not inside any module)
            modules.append(Module(a_domain=a_domain,
                                  domains_sequence=[DomainType.A]))
            a_domain_with_coords = next(a_domains_iter, None)
        elif a_domain_coords.start > module['components'][-1]['domain']['query_end']:  # module without A-domain
            modules.append(Module(a_domain=None,
                                  domains_sequence=get_domain_sequence(module)))
            module = next(modules_iter, None)
        else:  # A-domain inside the module
            modules.append(Module(a_domain=a_domain,
                                  domains_sequence=get_domain_sequence(module)))
            a_domain_with_coords = next(a_domains_iter, None)
            module = next(modules_iter, None)

    while module is not None:
        modules.append(Module(a_domain=None,
                              domains_sequence=get_domain_sequence(module)))
        module = next(modules_iter, None)
    while a_domain_with_coords is not None:
        a_domain, a_domain_coords = a_domain_with_coords
        modules.append(Module(a_domain=a_domain,
                              domains_sequence=[DomainType.A]))
        a_domain_with_coords = next(a_domains_iter, None)

    return modules


def check_orphan_c_domains(gene_data: dict,
                           config: antiSMASH_Processing_Config) -> Tuple[bool, bool]:  # (orphan_c_at_start, orphan_c_at_end)
    modules_start = min((domain['domain']['query_start']
                        for module in gene_data['modules']
                        for domain in module['components']),
                        default=None)

    modules_end = max((domain['domain']['query_end']
                      for module in gene_data['modules']
                      for domain in module['components']),
                        default=None)

    fst_c_domain_start = min((domain['query_start']
                             for domain in gene_data['domain_hmms']
                             if domain['hit_id'] in config.ANTISMASH_DOMAINS_NAMES_MAPPING and
                             DomainType[config.ANTISMASH_DOMAINS_NAMES_MAPPING[domain['hit_id']]].in_c_domain_group()),
                            default=None)
    last_c_domain_end = max((domain['query_end']
                            for domain in gene_data['domain_hmms']
                            if domain['hit_id'] in config.ANTISMASH_DOMAINS_NAMES_MAPPING and
                            DomainType[config.ANTISMASH_DOMAINS_NAMES_MAPPING[domain['hit_id']]].in_c_domain_group()),
                            default=None)

    orphan_c_start = fst_c_domain_start is not None and (
            modules_start is None or modules_start > fst_c_domain_start)
    orphan_c_end = last_c_domain_end is not None and (
            modules_end is None or modules_end < last_c_domain_end)
    return orphan_c_start, orphan_c_end


def get_a_domains_with_coords_per_gene(a_domains_info: Dict[A_Domain_Id, A_Domain],
                                       a_domains_coords: Dict[A_Domain_Id, Coords]) -> Dict[GeneId, List[Tuple[A_Domain, Coords]]]:
    a_domains_with_coords_per_gene = defaultdict(list)
    for a_domain_id, a_domain in sorted(a_domains_info.items(),
                                        key=lambda p: p[0].idx):  # sort by A-domain index so that they are added to the list in the order of appearance in the gene
        a_domains_with_coords_per_gene[a_domain_id.gene_id].append((a_domain, a_domains_coords[a_domain_id]))

    return a_domains_with_coords_per_gene


def extract_genes(contig_data: dict,
                  a_domains_info: Dict[A_Domain_Id, A_Domain],
                  a_domains_coords: Dict[A_Domain_Id, Coords],
                  config: antiSMASH_Processing_Config) -> List[Gene]:
    gene_coords = extract_gene_coords(contig_data)

    a_domains_with_coords_per_gene = get_a_domains_with_coords_per_gene(a_domains_info, a_domains_coords)

    genes = []
    gene_data_dict = contig_data['modules']['antismash.detection.nrps_pks_domains']['cds_results']
    for gene_id, gene_data in gene_data_dict.items():
        modules = extract_modules(gene_id, gene_data, a_domains_with_coords_per_gene[gene_id], config)
        orphan_c_at_start, orphan_c_at_end = check_orphan_c_domains(gene_data, config)
        genes.append(Gene(gene_id=gene_id,
                          coords=gene_coords[gene_id],
                          modules=modules,
                          orphan_c_at_start=orphan_c_at_start,
                          orphan_c_at_end=orphan_c_at_end))

    genes.sort(key=lambda gene: gene.coords.start)

    return list(filter(lambda gene: gene.modules, genes))


def extract_bgc_clusters(genome_id: str, ctg_idx: int,
                         contig_data: dict, genes: List[Gene]) -> List[BGC_Cluster]:
    bgcs = []
    for bgc_idx, bgc_data in enumerate(contig_data['areas'], start=1):
        if 'NRPS' in bgc_data['products']:
            bgc_genes = [gene for gene in genes
                         if bgc_data['start'] <= gene.coords.start <= gene.coords.end <= bgc_data['end']]

            if bgc_genes:
                bgcs.append(BGC_Cluster(bgc_id=BGC_ID(genome_id=genome_id,
                                                      contig_idx=ctg_idx,
                                                      bgc_idx=bgc_idx),
                                        genes=bgc_genes))
    return bgcs


def parse_antismash_json(antismash_json: antiSMASH_record,
                         config: antiSMASH_Processing_Config,
                         log: NerpaLogger) -> List[BGC_Cluster]:
    bgcs = []
    genome_id = antismash_json['input_file'].rsplit('.', 1)[0]  # remove extension
    for ctg_idx, contig_data in enumerate(antismash_json['records'], start=1):
        try:
            a_domains_info = extract_a_domains_info(contig_data, config,
                                                    ctg_idx, genome_id, log)  # pass genome_id and ctg_idx for error messages
            a_domains_coords = extract_a_domains_coords(contig_data)
            assert a_domains_info.keys() == a_domains_coords.keys()
            if not any(a_domain is not None for a_domain in a_domains_info.values()):  # no A-domains found in the contig
                continue

            genes = extract_genes(contig_data, a_domains_info, a_domains_coords, config)
            bgcs.extend(extract_bgc_clusters(genome_id, ctg_idx, contig_data, genes))
        except Exception as e:
            if config.DEBUG_MODE:
                log.error(f'Error parsing contig {ctg_idx} of genome {genome_id}:\n'
                          f'{type(e).__name__}: {e}')
            else:
                log.warning(f'Error parsing contig {ctg_idx} of genome {genome_id}')

    return bgcs
