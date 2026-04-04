# Project Overview — `nerpa_ms_core` (prototype v2)

## 0) One-liner
A Rust CLI tool that, for each `(HMM, rBAN_record, linearization)` item, builds a “nearby-variants” DAG over monomer sequences derived from a template molecule’s monomer graph (allowing limited edits). It then finds the top‑`N` most probable HMM state paths whose emitted monomer sequence matches some DAG path of total deviation weight `w ∈ {1,2,3}`, applies the implied edits to the template molecule, and outputs altered molecule variants as JSON.

## 1) Biological context (high-level)
We start from two complementary sources:

- **Genome-side prediction**: an **HMM** derived from genome analysis that encodes which peptide/monomer sequences are likely products of a biosynthetic gene cluster (BGC). The HMM assigns (log-)probability to candidate sequences, including non-emitting transitions where appropriate.
- **Molecule-side template**: a known or hypothesized peptide structure represented as a **Parsed_rBAN_Record**, i.e., a monomer graph describing the template molecule.

To relate these, we use an **NRP linearization**: an ordering of selected monomers in the template that provides an approximate correspondence between genome modules and positions along the template peptide.

The tool then **tweaks the template peptide** (via substitutions/insertions/deletions constrained by a deviation budget) to make it align better with the HMM, producing a ranked list of plausible “true product” variants. These variants are then tested downstream by another tool.

---

## 2) Inputs

Each input item is processed independently.

### 2.1 HMM
- Emitting + non-emitting states
- `state 0` = START, `state S-1` = FINISH
- Transition/emission scores are **log-probabilities** (`f64`)
- Non-emitting convention: `emissions[i].is_empty()`

```rust
#[derive(Debug, Clone, Deserialize)]
pub struct HMM {
    pub bgc_variant_id: BGC_Variant_ID,
    pub transitions: Vec<Vec<(StateIdx, LogProb)>>,
    pub emissions: Vec<Vec<LogProb>>,
}
```

The HMM emission alphabet corresponds to **MonomerCode** (used as DAG labels).

### 2.2 rBAN_Record (parsed template molecule)
```rust
#[derive(Debug, Clone, PartialEq)]
pub struct Parsed_rBAN_Record {
    pub compound_id: String,
    pub monomers: HashMap<MonomerIdx, MonomerInfo>,
    pub monomer_bonds: HashMap<MonomerEdge, MonomerEdgeInfo>,
    pub atoms: HashMap<AtomId, AtomInfo>,
    pub atomic_bonds: HashMap<AtomicEdge, AtomicEdgeInfo>,
    pub metadata: NRP_Metadata,
}
```

### 2.3 NRP_linearization
- `Vec<MonomerIdx>`
- Indices do not repeat and are a subset of those in the monomer graph.
- Intuition: a backbone-like ordering used to anchor correspondence between the template and the genome-side HMM, even if it conceptually joins multiple disjoint paths.

### 2.4 Global parameters
- `N`: number of top results to return per `w`
- weights: default `{1,2,3}`

---

## 3) Core internal types

### 3.1 DAG
A directed acyclic graph with:
- edges edges labeled with graph modifications
- vertices labeled with `Option<MonomerCode>` (the code of the monomer at that position, or `None` if it’s a non-emitting step)

A `START→FINISH` path in this DAG represents a candidate monomer sequence near the template/linearization, potentially involving edits.

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Edge<'a> {
    pub to: VertexId,
    pub weight: u8, // 0 or 1
    pub modification: Option<GraphModification<'a>>,
}

#[derive(Debug, Clone)]
pub struct VertexLabel {
	pub monomer_code: Option<MonomerCode>,
	pub name: String,  // used only for graph drawing
}

#[derive(Debug, Clone)]
pub struct DAG<'a> {
    pub nrp_variant_id: String,
    pub labels: Vec<VertexLabel>, // labels[v] = label of vertex v
    pub out_edges: Vec<Vec<Edge<'a>>>, // out_edges[v] = list of edges from vertex v
    pub start: VertexId,               // 0
    pub finish: VertexId,              // labels.len() - 1
}
```

### 3.2 Digitalized LogProb

To allow for a dynamic programming approach, I convert log probabilities to a discrete range -- integers from 0 to `MAX_DISCRETE_LOG_PROB` (e.g., 10000). The conversion is linear (except for rounding), mapping the smallest log-prob in the HMM to 0 and the largest to `MAX_DISCRETE_LOG_PROB`. I need an efficient bit-array like data structure with the following properties:
- Stores a set of `DISCRETE_LOG_PROB` values (integers in [0, `MAX_DISCRETE_LOG_PROB`])
- Supports bitwise shift and bitwise OR operations (for dynamic programming state updates)


---

## 4) DAG semantics: “nearby” sequences and deviation weights
A `START→FINISH` DAG path corresponds to a path in the underlying monomer graph with allowed deviations:
- substitutions
- insertions
- deletions

Edge weight meaning:
- **weight 0**: the DAG step mirrors an edge in the monomer graph (follows template connectivity)
- **weight 1**: the step represents a deviation from the monomer graph (edit/non-mirroring step)

Total path weight `w` is the sum of edge weights; the prototype computes results for `w ∈ {1,2,3}`.

---

## 5) Scoring / matching problem (HMM × DAG)

For each requested `w`:

1. Consider all DAG paths from START to FINISH with total deviation weight `w`.
2. Each DAG path spells a monomer sequence (the non-`None` vertex labels).
3. Find the **top‑`N` most probable HMM state paths** whose emitted sequence matches that monomer sequence.
   - HMM non-emitting states are epsilon steps (consume no DAG label).
4. Keep enough reconstruction info (internally) to later apply the solution back onto the template molecule.

---

## 6) Applying a solution back to the template molecule
Each selected solution implies surgical edits to the template monomer graph:

- substitute / insert / delete monomers inside the monomer graph (according to the chosen DAG path)
- produce a new `Parsed_rBAN_Record`
- produce an old→new monomer mapping for traceability

### Result struct (internal / output)
```rust
pub struct Altered_NRP_Variant {
    pub score: LogProb,
    pub new_molecule: Parsed_rBAN_Record,
    pub old_to_new_mon_map: Vec<(Option<MonomerIdx>, Option<MonomerIdx>)>,
}
```

---

## 7) Output (JSON, high-level)
For each input item:
- results grouped by deviation weight `w ∈ {1,2,3}`
- each group contains up to `N` `Altered_NRP_Variant` entries ranked by `score` (log-prob)

---

## 8) Tentative file structure (outdated, to be updated)

nerpa_ms_core/
├── Cargo.toml
├── README.md
├── fixtures/
│   ├── tiny_case.input.json
│   └── tiny_case.expected.json
├── src/
│   ├── main.rs
│   ├── lib.rs
│   ├── cli.rs
│   ├── types.rs
│   ├── io/
│   │   ├── mod.rs
│   │   ├── input.rs
│   │   ├── output.rs
│   │   └── json.rs
│   ├── hmm/
│   │   ├── mod.rs
│   │   ├── model.rs
│   │   └── validate.rs
│   ├── rban/
│   │   ├── mod.rs
│   │   ├── model.rs
│   │   ├── validate.rs
│   │   └── convert.rs
│   ├── monomer_graph/
│   │   ├── mod.rs
│   │   ├── model.rs
│   │   ├── build.rs
│   │   └── edit.rs
│   ├── dag/
│   │   ├── mod.rs
│   │   ├── model.rs
│   │   ├── build.rs
│   │   └── validate.rs
│   ├── solver/
│   │   ├── mod.rs
│   │   ├── candidate.rs
│   │   ├── topn.rs
│   │   └── viterbi_dag.rs
│   ├── apply/
│   │   ├── mod.rs
│   │   └── apply_candidate.rs
│   └── pipeline/
│       ├── mod.rs
│       └── run_item.rs
└── tests/
    ├── e2e_tiny_case.rs
    └── common/
        └── mod.rs


## 9) Plan and milestones

### Milestone 1 — Skeleton CLI + JSON input [DONE]
**Goal:** program builds, reads input JSON.

- CLI flags:
  - `--input PATH` (input JSON with the list of (HMM, Parsed_rBAN_Record, NRP_Linearization) triplets)
  - `--out PATH` (optional; else stdout)
  - `--max-edits` -- maximum number of edits allowed in a structure
  - `--num-variants-per-num-edits` -- the number of variants generated for each number of edits. 
- Implement input structs for the top-level JSON container.
- Parse and basic validation; print/emit a minimal per-item header in output.

**Deliverable:** `cargo run -- --input X.json --max-edits 3 --num-variants-per-num-edits 10` parses the input and prints a header for each item.

---

### Milestone 2 — Monomers database [DONE]
**Goal:** implement a database of monomers (from the rBAN monomer library) that can be used for surgical edits.

- Gather the rBAN monomers database: take parsed records from pnrpdb2 and take the most important monomers -- only supported residues + optional methylation.
- Identify a monomer's binding sites
- Function that, given a monomer, determines which other monomers it can be substituted for (based on binding site compatibility).

**Deliverable:** a monomer database module with a function `compatible_monomers(monomer) -> Vec<MonomerCode>` that can be used in the DAG construction step.

---

### Milestone 3 — DAG construction from (Parsed_rBAN_Record, linearization) [DONE]
**Goal:** build the deviation DAG with correct labels and 0/1 edge weights.

- Implement `build_dag(rban_record, linearization) -> DAG`.
- Ensure DAG vertex labels are `Option<MonomerCode>` and compatible with HMM’s `M`.
- Ensure DAG is acyclic (by construction or by check).

**Deliverable:** dump a debug summary (vertex/edge counts, maybe a few sample paths) and pass basic tests.

---

### Milestone 4 -- Discrete LogProb structure

**Goal:** implement the discrete log-probability structure and conversion function.

- Implement `DiscreteLogProb` data type with convertions to and from `LogProb`.
- Implement `DiscreteLogProbSet` (bit-array-like structure) with the operations: shift and OR operations.
   - `shift_left(&self, lp: LogProb) -> DiscreteLogProbSet` -- adds lp to all log probs values represented in the set. Values that go out of bounds are erased (like in usual bitwise shift of a bit-array)
   - `union(&self, other: DiscreteLogProbSet) -> DiscreteLogProbSet`


### Milestone 5 — HMM × DAG solver returning top‑N candidates (per weight)
**Goal:** compute ranked solutions *in the DAG/HMM world* (still no molecule edits).

- Implement `solve(hmm, dag, n, weights={1,2,3}) -> Vec<CandidatePerWeight>`.
- Candidate includes:
  - `score: LogProb`
  - backpointers for reconstruction (HMM states + DAG vertices + emitted MonomerCodes)
- Implement simple top‑N keeper (`Vec`, insert, sort, truncate).

**Deliverable:** on a tiny fixture, solver returns non-empty candidates for some `w`.

---

### Milestone 6 — Apply candidates to produce Altered_NRP_Variant
**Goal:** turn candidates into edited molecules + mapping.

- Implement `apply_candidate(rban_record, candidate) -> Altered_NRP_Variant`.
- Keep the edits “surgical” and localized (prototype scope).
- Produce `old_to_new_mon_map`.

**Deliverable:** end-to-end run produces actual `new_molecule` records for returned candidates.

---

### Milestone 7 — Tests and fixtures
**Goal:** keep regressions from sneaking in.

- Unit tests:
  - HMM validation
  - linearization validation
  - DAG build invariants (start/finish labels, edge weights in {0,1}, label range)
- End-to-end test on a tiny synthetic case:
  - 1 item, small HMM, small monomer graph, small linearization
  - checks: output groups exist, ranking order by score, stable serialization shape

**Deliverable:** `cargo test` runs meaningful coverage for parsing + DAG build + solver + apply.

---

### Milestone 8 (optional) — Debuggability and performance niceties
**Goal:** make it easier to trust and iterate.

- Optional `--debug` output fields (hmm_states, dag_vertices, emitted codes).
- Better top‑N structure (binary heap) if needed.
- Timing logs per stage (parsing, DAG build, solve, apply).

**Deliverable:** easier diagnosis on real inputs without changing core logic.

