use crate::data_types::{dag::{DAG, Edge}, discrete_log_prob::DiscreteLogProb, dp_table::{BacktrackPointer, DP_Coords, DP_Table}, hmm::StateIdx};

pub struct Solution <'a>{
    pub states: Vec<StateIdx>,
    pub dag_edges: Vec<Edge<'a>>,
    pub dlp: DiscreteLogProb,
    pub weight: usize,
}


fn is_dp_start(coords: &DP_Coords) -> bool {
    coords.vertex == 0 && coords.weight == 0 && coords.state == 0
}

#[derive(Debug, Clone)]
struct Frame<'a> {
    coords: DP_Coords,
    dlp: DiscreteLogProb,
    // built backwards (finish -> start)
    states_rev: Vec<StateIdx>,
    dag_edges_rev: Vec<Edge<'a>>,
}

impl<'a> Frame<'a> {
    fn into_solution(mut self, weight: usize) -> Solution<'a> {
        self.states_rev.reverse();
        self.dag_edges_rev.reverse();

        Solution {
            states: self.states_rev,
            dag_edges: self.dag_edges_rev,
            dlp: self.dlp,
            weight,
        }
    }
}

pub struct BacktrackSolutionsIter<'a> {
    dp: &'a DP_Table<'a>,
    weight: usize,
    end_coords: DP_Coords,

    // Discrete log-probs available at (DAG_FINISH, weight, HMM_FINISH), descending.
    end_dlp_iter: Box<dyn Iterator<Item = DiscreteLogProb> + 'a>,

    // DFS stack for the current end_dlp
    stack: Vec<Frame<'a>>,
}

impl<'a> BacktrackSolutionsIter<'a> {
    fn start_new_dlp(&mut self, end_dlp: DiscreteLogProb) {
        self.stack.clear();
        self.stack.push(Frame {
            coords: self.end_coords,
            dlp: end_dlp,
            states_rev: vec![self.end_coords.state],
            dag_edges_rev: Vec::new(),
        });
    }

    fn expand_one(&mut self, frame: Frame<'a>) {
        for (parent, parent_dlp, dag_edge) in self.dp.get_backtrack_parents(&frame.coords, frame.dlp) {
            let mut next = frame.clone();
            next.coords = parent;
            next.dlp = parent_dlp;

            // we moved to the parent -> record parent's state
            next.states_rev.push(parent.state);

            // record DAG edge if this DP transition corresponds to a DAG step
            if let Some(e) = dag_edge {
                next.dag_edges_rev.push(e);
            }

            self.stack.push(next);
        }
    }
}

impl<'a> Iterator for BacktrackSolutionsIter<'a> {
    type Item = Solution<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some(frame) = self.stack.pop() {
                if is_dp_start(&frame.coords) {
                    return Some(frame.into_solution(self.weight));
                }
                self.expand_one(frame);
                continue;
            }

            let next_end_dlp = self.end_dlp_iter.next()?;
            self.start_new_dlp(next_end_dlp);
        }
    }
}

pub fn backtrack_solutions<'a>(
    weight: usize,
    dp: &'a DP_Table<'a>,
    dag: &'a DAG<'a>,
) -> BacktrackSolutionsIter<'a> {
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
        stack: Vec::new(),
    }
}
