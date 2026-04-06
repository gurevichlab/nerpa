use super::{common_types::LogProb, discrete_log_prob::DiscreteLogProbSet};


// A small wrapper around a flat Vec
// Logical layout: dp[vertex][weight][state]
#[derive(Debug, Clone)]
pub struct DP_Table {
    n_vertices: usize,
    n_weights: usize, // = max_weight + 1
    n_states: usize,
    data: Vec<DiscreteLogProbSet>,
    parents: Vec<Vec<(usize, Option<LogProb>)>>, // parallel to data, stores parent indices and optional shift for each cell
}

impl DP_Table {
    pub fn new(n_vertices: usize, max_weight: usize, n_states: usize) -> Self {
        let n_weights = max_weight + 1;
        let len = n_vertices * n_weights * n_states;

        let mut data = vec![DiscreteLogProbSet::empty(); len];
	// start state has probability 1
	data[0] = DiscreteLogProbSet::from_logprob_vec(vec![0.0]);

	let parents = vec![Vec::new(); len];

        Self {
            n_vertices,
            n_weights,
            n_states,
            data,
	    parents
        }
    }

    #[inline]
    pub fn idx(&self, vertex: usize, weight: usize, state: usize) -> usize {
        debug_assert!(vertex < self.n_vertices);
        debug_assert!(weight < self.n_weights);
        debug_assert!(state < self.n_states);

        (vertex * self.n_weights + weight) * self.n_states + state
    }

    #[inline]
    pub fn get(&self, vertex: usize, weight: usize, state: usize) -> &DiscreteLogProbSet {
        let i = self.idx(vertex, weight, state);
        &self.data[i]
    }

    #[inline]
    pub fn get_parents(&self, vertex: usize, weight: usize, state: usize) -> &Vec<(usize, Option<LogProb>)> {
        let i = self.idx(vertex, weight, state);
        &self.parents[i]
    }

    pub fn n_vertices(&self) -> usize {
        self.n_vertices
    }

    pub fn max_weight(&self) -> usize {
        self.n_weights - 1
    }

    pub fn n_states(&self) -> usize {
        self.n_states
    }

    #[inline]
    fn update_cell_with_cell(&mut self, dst: usize, src: usize) {
        if dst == src {
            return; // union with itself does nothing
        }

        // Split so we can borrow two different elements safely.
        if src < dst {
            let (left, right) = self.data.split_at_mut(dst);
            let src_ref = &left[src];
            let dst_ref = &mut right[0];
            dst_ref.union_inplace(src_ref);
        } else {
            let (left, right) = self.data.split_at_mut(src);
            let dst_ref = &mut left[dst];
            let src_ref = &right[0];
            dst_ref.union_inplace(src_ref);
        }
    }

    fn update_from_idx(&mut self, dst: usize, src: usize, shift: Option<LogProb>) {
	match shift {
	    Some(s) => {
		// Shift src by s and union into dst.
		let shifted_src = self.data[src].add_to_all(s);
		self.data[dst].union_inplace(&shifted_src);
	    },
	    None => {
		// No shift, just union src into dst.
		self.update_cell_with_cell(dst, src);
	    }
	}

	self.parents[dst].push((src, shift));
    }

    // Update dst cell by unioning in src cell, optionally applying a shift to all log-probabilities from src before unioning.
    pub fn update(&mut self,
		  dst: (usize, usize, usize),
		  src: (usize, usize, usize),
		  shift: Option<LogProb>) {
	let src_idx = self.idx(src.0, src.1, src.2);
	let dst_idx = self.idx(dst.0, dst.1, dst.2);
	
	self.update_from_idx(src_idx, dst_idx, shift);
    }
}

