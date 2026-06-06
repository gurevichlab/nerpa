use crate::data_types::{mod_graph::{ModGraph, Edge, VertexId}, discrete_log_prob::DiscreteLogOdds, dp_table::{BacktrackPointer, DP_Coords, DP_Table}, hmm::StateIdx};

pub struct Solution <'mon_db>{
    pub states: Vec<StateIdx>,
    pub dag_edges: Vec<Edge<'mon_db>>,
    pub dlo: DiscreteLogOdds,
    pub weight: usize,
}

fn is_dp_start(coords: &DP_Coords) -> bool {
    coords.vertex == 0 && coords.weight == 0 && coords.state == 0
}

#[derive(Debug, Clone)]
struct Frame<'mon_db> {
    coords: DP_Coords,
    dlo: DiscreteLogOdds,
    // built backwards (finish -> start)
    states_rev: Vec<StateIdx>,
    dag_edges_rev: Vec<Edge<'mon_db>>,
}

impl<'mon_db> Frame<'mon_db> {
    fn into_solution(mut self, weight: usize, dlo: DiscreteLogOdds) -> Solution<'mon_db> {
        self.states_rev.reverse();
        self.dag_edges_rev.reverse();

        Solution {
            states: self.states_rev,
            dag_edges: self.dag_edges_rev,
            dlo,
            weight,
        }
    }
}

pub struct BacktrackSolutionsIter<'mon_db: 'iter, 'iter> {
    dp: &'iter DP_Table<'mon_db>,
    weight: usize,
    end_coords: DP_Coords,

    // Discrete log-probs available at (DAG_FINISH, weight, HMM_FINISH), descending.
    end_dlo_iter: Box<dyn Iterator<Item = DiscreteLogOdds> + 'iter>,
    cur_dlo: DiscreteLogOdds,

    // DFS stack for the current end_dlp
    stack: Vec<Frame<'mon_db>>,
}

use crate::data_types::discrete_log_prob::SCALING_FACTOR;
use crate::algo::generic::rounded;

impl<'mon_db, 'iter> BacktrackSolutionsIter<'mon_db, 'iter> {
    fn start_new_dlp(&mut self, end_dlo: DiscreteLogOdds) {
        self.stack.clear();
        self.stack.push(Frame {
            coords: self.end_coords,
            dlo: end_dlo,
            states_rev: vec![self.end_coords.state],
            dag_edges_rev: Vec::new(),
        });
	self.cur_dlo = end_dlo;
    }

    fn expand_one(&mut self, frame: Frame<'mon_db>) {
	let debug = false;
	if debug {
	    let lp_rounded = rounded(frame.dlo.to_logodds(), 2);
	    println!("Expanding frame:\n\tcoords={:?}\n\tdlp={},\tlp={}\n\tstates_rev={:?}\n\tdag_edges_rev={:?}", frame.coords, frame.dlo.0, lp_rounded, frame.states_rev, frame.dag_edges_rev);
	    println!("Number of backtracking pointers: {}", self.dp.get_backtrack_pointers(&frame.coords).len());
	}

        for ptr in self.dp.get_backtrack_pointers(&frame.coords) {
	    let parent_dlo = match ptr.dlo_shift {
		// shift is applied parent -> child
		// so we need to apply -shift to get from child back to parent
		Some(shift) => {
		    if shift == isize::MIN {
			None  // avoid overflow when negating shift
		    }
		    else {
			frame.dlo.shift(-shift)
		    }
		},
		None => Some(frame.dlo),
	    };
	    if parent_dlo.is_none() || !self.dp.get(&ptr.parent).contains(parent_dlo.unwrap()) { continue }

	    if debug {
		let lo_shift = ptr.dlo_shift
		    .map(|lp| rounded((lp as f64) / SCALING_FACTOR, 2));
		println!("Backtracking pointer:\n\tparent={:?}\n\tdlp_shift={:?}\tlp={:?}\n\tdag_edge={:?}", ptr.parent, ptr.dlo_shift, lo_shift, ptr.dag_edge);
	    }
            let mut next = frame.clone();
            next.coords = ptr.parent;
            next.dlo = parent_dlo.unwrap();

            // DP transition corresponds to a change of HMM state
	    // -> append parent's state
	    if ptr.dlo_shift.is_some() {
		next.states_rev.push(ptr.parent.state);
	    }

            // DP transition corresponds to a DAG step
	    // append the corresponding DAG edge
            if let Some(e) = ptr.dag_edge {
                next.dag_edges_rev.push(e);
            }

            self.stack.push(next);
        }
    }
}

impl<'mon_db, 'iter> Iterator for BacktrackSolutionsIter<'mon_db, 'iter> {
    type Item = Solution<'mon_db>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(frame) = self.stack.pop() {
		//println!("Backtracking frame: coords={:?}, states_rev={:?}, dag_edges_rev={:?}", frame.coords, frame.states_rev, frame.dag_edges_rev);
                if is_dp_start(&frame.coords) {
                    return Some(frame.into_solution(self.weight, self.cur_dlo));
                }
                self.expand_one(frame);
            }
	    else {
		let next_end_dlp = self.end_dlo_iter.next()?;
		self.start_new_dlp(next_end_dlp);
	    }
        }
    }
}

pub fn backtrack_solutions<'mon_db: 'iter, 'iter>(
    weight: usize,
    dp: &'iter DP_Table<'mon_db>,
    dag: &'iter ModGraph<'mon_db>,
) -> BacktrackSolutionsIter<'mon_db, 'iter> {
    debug_assert!(weight <= dp.max_weight());

    let end_coords = DP_Coords {
        vertex: dag.finish,
        weight,
        state: dp.n_states() - 1, // HMM FINISH
    };

    let end_dlo_set = dp.get(&end_coords);

    BacktrackSolutionsIter {
        dp,
        weight,
        end_coords,
        end_dlo_iter: Box::new(end_dlo_set.iter_desc()),
	cur_dlo: DiscreteLogOdds::from_logodds(0.0),  // doesn't matter
        stack: Vec::new(),
    }
}
