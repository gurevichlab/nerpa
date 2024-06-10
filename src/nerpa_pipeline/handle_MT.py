#!/usr/bin/env python

import src.nerpa_pipeline.handle_helper as handle_helper

def get_modules_idxs_with_methylation(gene: Gene) -> List[int]:
    return []


def get_MT_AA(dirname, orfs_order):
    MT_AA_list = []
    # ctg_orf, ctg_orf_Aid, domain_type(AMP-binding, PCP, MT)
    domains=handle_helper.get_domains_list(dirname)
    orf_ori = handle_helper.get_orf_orientation(dirname)

    # reverse domains on - strand
    for i in range(len(domains)):
        if domains[i]:
            if orf_ori[domains[i][0][0]] == '-':
                domains[i].reverse()

    order_dom = []
    for corf in orfs_order:
        for j in range(len(domains)):
            if domains[j][0][0] == corf:
                order_dom.append(domains[j])
    domains = order_dom
    is_AA = lambda dlst : ("PKS" in dlst[2]) or ("AMP-binding" == dlst[2])

    for orfds in domains:
        cur_i = 0
        has_mt = False
        for i in range(len(orfds)):
            if is_AA(orfds[i]):
                if has_mt:
                    MT_AA_list.append(orfds[cur_i][0] + "_" + orfds[cur_i][1].split('_')[-1])

                cur_i = i
                has_mt = False

            if orfds[i][-1] == "MT":
                has_mt = True

        if has_mt:
            MT_AA_list.append(orfds[cur_i][0] + "_" + orfds[cur_i][1].split('_')[-1])

    return MT_AA_list
