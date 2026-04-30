use crate::data_types::{dag::{DAG, Edge, VertexId}, discrete_log_prob::DiscreteLogProb, dp_table::{BacktrackPointer, DP_Coords, DP_Table}, hmm::StateIdx};

pub struct Solution <'mon_db>{
    pub states: Vec<StateIdx>,
    pub dag_edges: Vec<Edge<'mon_db>>,
    pub dlp: DiscreteLogProb,
    pub weight: usize,
}

fn is_dp_start(coords: &DP_Coords) -> bool {
    coords.vertex == 0 && coords.weight == 0 && coords.state == 0
}

#[derive(Debug, Clone)]
struct Frame<'mon_db> {
    coords: DP_Coords,
    dlp: DiscreteLogProb,
    // built backwards (finish -> start)
    states_rev: Vec<StateIdx>,
    dag_edges_rev: Vec<Edge<'mon_db>>,
}

impl<'mon_db> Frame<'mon_db> {
    fn into_solution(mut self, weight: usize, dlp: DiscreteLogProb) -> Solution<'mon_db> {
        self.states_rev.reverse();
        self.dag_edges_rev.reverse();

        Solution {
            states: self.states_rev,
            dag_edges: self.dag_edges_rev,
            dlp,
            weight,
        }
    }
}

pub struct BacktrackSolutionsIter<'mon_db: 'iter, 'iter> {
    dp: &'iter DP_Table<'mon_db>,
    weight: usize,
    end_coords: DP_Coords,

    // Discrete log-probs available at (DAG_FINISH, weight, HMM_FINISH), descending.
    end_dlp_iter: Box<dyn Iterator<Item = DiscreteLogProb> + 'iter>,
    cur_dlp: DiscreteLogProb,

    // DFS stack for the current end_dlp
    stack: Vec<Frame<'mon_db>>,
}

use crate::data_types::discrete_log_prob::SCALING_FACTOR;

pub fn rounded(f: f64, digits: usize) -> f64 {
    if !f.is_finite() {
	return f;
    }
    if digits == 0 {
	return f.round();
    }
    let factor = 10f64.powf(digits as f64);
    if !factor.is_finite() || factor == 0.0 {
	return f;
    }
    (f * factor).round() / factor
}

impl<'mon_db, 'iter> BacktrackSolutionsIter<'mon_db, 'iter> {
    fn start_new_dlp(&mut self, end_dlp: DiscreteLogProb) {
        self.stack.clear();
        self.stack.push(Frame {
            coords: self.end_coords,
            dlp: end_dlp,
            states_rev: vec![self.end_coords.state],
            dag_edges_rev: Vec::new(),
        });
	self.cur_dlp = end_dlp;
    }

    fn expand_one(&mut self, frame: Frame<'mon_db>) {
	let debug = false;
	if debug {
	    let lp_rounded = rounded(frame.dlp.to_logprob(), 2);
	    println!("Expanding frame:\n\tcoords={:?}\n\tdlp={},\tlp={}\n\tstates_rev={:?}\n\tdag_edges_rev={:?}", frame.coords, frame.dlp.0, lp_rounded, frame.states_rev, frame.dag_edges_rev);
	    println!("Number of backtracking pointers: {}", self.dp.get_backtrack_pointers(&frame.coords).len());
	}

	let ptrs = self.dp
	    .get_backtrack_pointers(&frame.coords) 
	    .iter()
	    .filter(|ptr| {
		let parent_dlp = match ptr.dlp_shift {
		    Some(shift) => frame.dlp.shift(shift as i64),
		    None => Some(frame.dlp),
		};
		if let Some(dlp) = parent_dlp {
		    self.dp.get(&ptr.parent).contains(dlp)
		}
		else { false }
	    });
        for ptr in self.dp.get_backtrack_pointers(&frame.coords) {
	    let parent_dlp = match ptr.dlp_shift {
		Some(shift) => frame.dlp.shift(shift as i64),
		None => Some(frame.dlp),
	    };
	    if parent_dlp.is_none() || !self.dp.get(&ptr.parent).contains(parent_dlp.unwrap()) { continue }

	    if debug {
		let lp_shift = ptr.dlp_shift
		    .map(|lp| rounded((lp as f64) / SCALING_FACTOR, 2));
		println!("Backtracking pointer:\n\tparent={:?}\n\tdlp_shift={:?}\tlp={:?}\n\tdag_edge={:?}", ptr.parent, ptr.dlp_shift, lp_shift, ptr.dag_edge);
	    }
            let mut next = frame.clone();
            next.coords = ptr.parent;
            next.dlp = parent_dlp.unwrap();

            // we moved to the parent -> record parent's state
	    if ptr.dlp_shift.is_some() {
		next.states_rev.push(ptr.parent.state);
	    }

            // record DAG edge if this DP transition corresponds to a DAG step
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
                    return Some(frame.into_solution(self.weight, self.cur_dlp));
                }
                self.expand_one(frame);
            }
	    else {
		let next_end_dlp = self.end_dlp_iter.next()?;
		self.start_new_dlp(next_end_dlp);
	    }
        }
    }
}

pub fn backtrack_solutions<'mon_db: 'iter, 'iter>(
    weight: usize,
    dp: &'iter DP_Table<'mon_db>,
    dag: &'iter DAG<'mon_db>,
) -> BacktrackSolutionsIter<'mon_db, 'iter> {
    debug_assert!(weight <= dp.max_weight());

    let end_coords = DP_Coords {
        vertex: dag.finish,
        weight,
        state: dp.n_states() - 1, // HMM FINISH
    };

    let end_dlp_set = dp.get(&end_coords);

    BacktrackSolutionsIter {
        dp,
        weight,
        end_coords,
        end_dlp_iter: Box::new(end_dlp_set.iter_desc()),
	cur_dlp: DiscreteLogProb::from_logprob(0.0),  // doesn't matter
        stack: Vec::new(),
    }
}
