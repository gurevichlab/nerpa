use crate::data_types::dag::{DAG, Edge};
use crate::data_types::graph_modifications::{get_possible_subs, GraphModification};
use crate::data_types::{monomers_db::MonomersDB, parsed_rban_record::Parsed_rBAN_Record};

pub fn create_dag<'a>(rban_record: &Parsed_rBAN_Record, monomers_db: &'a MonomersDB) -> DAG<'a> {
    let mut labels: Vec<Option<String>> = vec![Some("INITIAL".to_string())];
    let mut out_edges: Vec<Vec<Edge>> = vec![Vec::new()];
    let mut subgraph_root = 0;

    for (monomer_idx, monomer_info) in rban_record.monomers.iter() {
        let subs = {
            let subs = get_possible_subs(rban_record, (*monomer_idx).clone(), monomers_db);
            // filter out substitutions that are identical to the original monomer
            subs.into_iter()
                .take(3) // for debugging purposes, limit the number of substitutions to 3
                .filter(|gm| {
                    if let GraphModification::Substitute {
                        monomer_idx: _,
                        mon_db_entry,
                    } = gm
                    {
                        mon_db_entry.monomer.name != monomer_info.name
                    } else {
                        true
                    }
                })
                .collect::<Vec<_>>()
        };
        let next_subgraph_root = subgraph_root + subs.len() + 1;

        labels.push(Some(monomer_info.name.0.clone() + "*"));
        out_edges[subgraph_root].push(Edge {
            to: subgraph_root + 1,
            weight: 0,
            modification: None,
        });

        assert_eq!(out_edges.len(), subgraph_root + 1);
        out_edges.push(vec![Edge {
            to: next_subgraph_root,
            weight: 0,
            modification: None,
        }]);

        for (i, gm) in subs.into_iter().enumerate() {
            match &gm {
                GraphModification::Substitute { mon_db_entry, .. } => {
                    labels.push(Some(mon_db_entry.monomer.name.0.clone()));
                }
                _ => unreachable!("expected Substitute"),
            }

            out_edges[subgraph_root].push(Edge {
                to: subgraph_root + 1 + i,
                weight: 1,
                modification: Some(gm),
            });
            out_edges.push(vec![Edge {
                to: next_subgraph_root,
                weight: 0,
                modification: None,
            }]);
        }

        subgraph_root = next_subgraph_root;
	labels.push(None);
	out_edges.push(Vec::new());
    }

    let labels_len = labels.len();
    labels[labels_len - 1] = Some("FINAL".to_string());
    DAG { 
	nrp_variant_id: rban_record.compound_id.clone(),
	labels,
	out_edges,
	start: 0,
	finish: labels_len - 1,
    }
}
