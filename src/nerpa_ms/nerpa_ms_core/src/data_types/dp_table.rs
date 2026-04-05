use super::discrete_log_prob::DiscreteLogProbSet;


// A small wrapper around a flat Vec
// Logical layout: dp[vertex][weight][state]
#[derive(Debug, Clone)]
pub struct DP_Table {
    n_vertices: usize,
    n_weights: usize, // = max_weight + 1
    n_states: usize,
    data: Vec<DiscreteLogProbSet>,
}

impl DP_Table {
    pub fn new(n_vertices: usize, max_weight: usize, n_states: usize) -> Self {
        let n_weights = max_weight + 1;
        let len = n_vertices * n_weights * n_states;

        // Assumes DiscreteLogProbSet::empty() exists.
        // If your type uses Default instead, replace with `DiscreteLogProbSet::default()`.
        let data = vec![DiscreteLogProbSet::empty(); len];

        Self {
            n_vertices,
            n_weights,
            n_states,
            data,
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
    pub fn get_mut(
        &mut self,
        vertex: usize,
        weight: usize,
        state: usize,
    ) -> &mut DiscreteLogProbSet {
        let i = self.idx(vertex, weight, state);
        &mut self.data[i]
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
    pub fn union_cell_into_cell(&mut self, src: usize, dst: usize) {
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
}

