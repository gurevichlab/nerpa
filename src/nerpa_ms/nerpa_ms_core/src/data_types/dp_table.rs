use super::{common_types::LogProb, dag::Edge, discrete_log_prob::{DiscreteLogProb, DiscreteLogProbSet}};

pub type DP_Idx = usize;

#[derive(Debug, Clone, Copy)]
pub struct BacktrackPointer<'a>{
    parent: DP_Coords,
    dlp_shift: usize, // shift to apply to all log-probs from parent cell when backtracking
    dag_edge: Option<Edge<'a>>,
}

// A small wrapper around a flat Vec
// Logical layout: dp[vertex][weight][state]
#[derive(Debug, Clone)]
pub struct DP_Table <'a>{
    n_vertices: usize,
    n_weights: usize, // = max_weight + 1
    n_states: usize,
    data: Vec<DiscreteLogProbSet>,
    parents: Vec<Vec<BacktrackPointer<'a>>>, // parallel to data, stores parent indices and optional shift for each cell
}

#[derive(Debug, Clone, Copy)]
pub struct DP_Coords {
    pub vertex: usize,
    pub weight: usize,
    pub state: usize,
}

impl <'a> DP_Table<'a>{
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
    pub fn idx(&self, coords: &DP_Coords) -> DP_Idx {
	let &DP_Coords { vertex, weight, state } = coords;
        debug_assert!(vertex < self.n_vertices);
        debug_assert!(weight < self.n_weights);
        debug_assert!(state < self.n_states);

        (vertex * self.n_weights + weight) * self.n_states + state
    }

    #[inline]
    pub fn idx_to_coordinates(&self, idx: DP_Idx) -> DP_Coords {
	debug_assert!(idx < self.data.len());

	let vertex = idx / (self.n_weights * self.n_states);
	let weight = (idx / self.n_states) % self.n_weights;
	let state = idx % self.n_states;
	DP_Coords{vertex, weight, state}
    }

    #[inline]
    pub fn get(&self, coords: &DP_Coords) -> &DiscreteLogProbSet {
        let i = self.idx(coords);
        &self.data[i]
    }

    #[inline]
    pub fn get_by_idx(&self, idx: DP_Idx) -> &DiscreteLogProbSet {
	debug_assert!(idx < self.data.len());
        &self.data[idx]
    }

    #[inline]
    pub fn get_backtrack_parents(&self, coords: &DP_Coords, dlp: DiscreteLogProb) -> Vec<(DP_Coords, DiscreteLogProb, Option<Edge<'a>>)> {
	let idx = self.idx(&coords);
	debug_assert!(idx < self.data.len());
	debug_assert!(self.data[idx].contains(dlp), "get_backtrack_parents: requested dlp not in cell");

	self.parents[idx]
	    .iter()
	    .filter_map(|&ptr| {
		let BacktrackPointer{parent: parent,
				     dlp_shift: dlp_shift,
				     dag_edge: dag_edge} = ptr;
		let parent_idx = self.idx(&parent);
		let parent_dlp = dlp.shift(dlp_shift as i64);
		if self.data[parent_idx].contains(parent_dlp) {
		    Some((parent, parent_dlp, dag_edge))
		} else {None}
	    })
	    .collect()
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

    fn update_from_idx(&mut self, dst: usize, src: usize, shift: Option<LogProb>, dag_edge: Option<Edge<'a>>) {
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

	let dlp_shift = if let Some(lp) = shift {
	    DiscreteLogProbSet::get_abs_shift(lp)
	} else {0};

	let parent = self.idx_to_coordinates(src);
	self.parents[dst].push(BacktrackPointer{
	    parent,
	    dlp_shift,
	    dag_edge,
	});
    }

    // Update dst cell by unioning in src cell, optionally applying a shift to all log-probabilities from src before unioning.
    pub fn update(&mut self,
		  dst: &DP_Coords,
		  src: &DP_Coords,
		  shift: Option<LogProb>,
		  dag_edge: Option<Edge<'a>>) {
	let src_idx = self.idx(src);
	let dst_idx = self.idx(dst);
	
	self.update_from_idx(dst_idx, src_idx, shift, dag_edge);
    }
}

