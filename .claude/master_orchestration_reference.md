<!-- REFERENCE COPY — DO NOT EDIT DIRECTLY
     Source of truth: trasim-lab/MASTER_CLAUDE.md
     Active copy loaded via: .claude/rules/master_orchestration.md (symlink, gitignored)
     This file is git-tracked for safekeeping. To update, re-copy from the source.

# MASTER_CLAUDE.md – Botos Transcription Simulation Project

This file contains the **shared coding rules** and **master orchestration protocol** for the Botos Transcription Simulation Project. All Claude agents operating in any tool repository (POLYARIS, KINEMA, GENEFIT, SCAN, ProtoPos) **must** read this file alongside their tool-specific CLAUDE.md.

**Location:** `trasim-lab/MASTER_CLAUDE.md`

---

## Architecture

trasim-lab is the orchestration layer for the project. It coordinates multiple Claude Code instances working across five tool repositories.

**Dependency direction (inviolable):**

```
trasim-lab → POLYARIS
trasim-lab → KINEMA
trasim-lab → GENEFIT
trasim-lab → SCAN
trasim-lab → ProtoPos
```

**Never:** A tool repo must NEVER import or depend on trasim-lab. That would destroy modularity. Lab orchestration lives exclusively in trasim-lab.

**Directory structure:**

```
trasim-lab/
├── MASTER_CLAUDE.md              # This file (master orchestration rules)
├── README.md                     # Architecture overview
├── pyproject.toml                # Tool repos as dependencies
├── agents/                       # Agent governance and coordination
│   ├── multi_agent_work_plan.md  # Full governance specification (v3)
│   ├── agent_lab_operator_guide.docx
│   ├── generate_operator_guide.py
│   ├── status_board.md           # Current task tracker
│   ├── research_suggestions.md   # End-of-session suggestions
│   ├── checkpoints/              # Session checkpoint reports
│   ├── handoffs/                 # Cross-repo handoff notes
│   ├── reviews/                  # Reviewer decisions
│   └── collector_output/         # Aggregated reports for review
├── workflows/                    # Cross-repo automation scripts
│   ├── collect_reports.sh        # Two-pass reviewer driver
│   └── nightly_smoke_test.sh     # Integration smoke test
├── configs/                      # Cross-repo specifications
│   ├── interface_schema.yaml     # Machine-level API contracts (semver)
│   └── living_spec.md            # Architectural intent document
├── evaluation/                   # Research tracking
│   ├── research_goals.md
│   ├── research_plan.md
│   ├── framework_description.md
│   └── interesting_results/      # Cross-repo findings
└── logs/                         # Automation logs
    └── smoke_test_*.log
```

---

## Python Version Support

The project targets **Python 3.10** as the primary development version. Python **3.10, 3.11, 3.12, and 3.13** are all supported.

---

## Conda Environment Setup

### Unified Environment (Recommended)

A single unified conda environment covers **all five** packages. This is the recommended setup for most users.

```bash
# Quick setup (from the project root):
bash setup_environment.sh              # Default: Python 3.10
bash setup_environment.sh 3.12         # Or specify another version: 3.11, 3.12, 3.13

# Manual setup:
conda env create -f environment.yml
conda activate polyaris_310
pip install -e POLYARIS -e KINEMA -e GENEFIT -e SCAN
# ProtoPos is imported directly (no pip install needed).
```

### Per-Repository Environments (Alternative)

Each repository also has its own dedicated conda environment pinned to Python 3.10. These are useful when working on a single repo in isolation.

| Repository | Conda Environment | Environment File |
|------------|-------------------|------------------|
| KINEMA     | `kinema_310`      | `KINEMA/environment.yml` |
| POLYARIS   | `polyaris_310`    | `POLYARIS/environment.yml` |
| GENEFIT    | `genefit_310`     | `GENEFIT/environment.yml` |
| SCAN       | `scan_310`        | `SCAN/environment.yml` |
| ProtoPos   | `protopos`        | *(standalone conda env)* |

### Conda Environment Activation

Before running **any** command (tests, scripts, pipelines, model generation, etc.), Claude **must** activate a conda environment first.

**Preferred:** Use the unified environment for all work:

```bash
source ~/miniconda3/bin/activate polyaris_310
```

**Alternative:** Use a per-repo environment:

```bash
source ~/miniconda3/bin/activate <env_name>
```

**Rules:**
- Always activate the correct environment before running Python, pytest, pip, or any other command that depends on the repo's packages.
- When working across multiple repos in a single shell command, use the unified environment.
- Never assume the correct environment is already active. Always include the activation in your command.

---

## Script Path Resolution Policy

All executable scripts (shell scripts and Python entry points) **must** work correctly when invoked from **any** working directory, not just from the directory where the script lives.

**Required patterns:**

For **shell scripts** (`.sh`):
```bash
# At the top of every script, resolve absolute paths from the script's own location.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$script_dir/.." && pwd)"

# Then use $REPO_DIR or $script_dir to build all paths.
# GOOD:
RESULTS_DIR="${REPO_DIR}/results"
CONFIG_DIR="${REPO_DIR}/config_files"

# BAD (breaks when script is called from a different directory):
RESULTS_DIR="../results"
CONFIG_DIR="../config_files"
```

For **Python files** (`.py`):
```python
# Resolve the repository root relative to __file__, not os.getcwd().
# GOOD:
repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
logs_dir = os.path.join(repo_root, "logs")

# BAD (breaks when script is called from a different directory):
logs_dir = os.path.join("..", "logs")
```

**Rules:**
- **Never** use bare relative paths like `../results`, `../logs`, `../config_files`, or `../screenshots` in executable code.
- **Always** derive paths from the script's own location using `${BASH_SOURCE[0]}` (shell) or `__file__` (Python).
- Relative paths in comments, docstrings, help text examples, and non-path string literals are fine.
- This rule applies to all repositories: POLYARIS, KINEMA, GENEFIT, SCAN, and ProtoPos.

---

## Bug Fixing Policy

When debugging an error, Claude **must** find and fix the **root cause**. Claude must **never** make the bug silent by adding graceful handling, try/except wrappers, default returns, or guard clauses that suppress the symptom without understanding why it occurs.

**Prohibited patterns:**
- Adding `if len(x) == 0: return` to hide a mismatch that should not happen.
- Wrapping a crash site in `try/except` to swallow the error.
- Adding fallback values or empty-array guards at the symptom location instead of tracing why the data is wrong upstream.
- Any modification to the function that **receives** bad data, when the real fix is in the code that **produces** bad data.

**Required approach:**
1. When an error occurs, trace the data flow **upstream** to find where the incorrect state originates.
2. Add temporary debug prints or assertions to confirm the root cause before applying a fix.
3. Fix the code that **produces** the wrong data, not the code that **consumes** it.
4. Remove all debug prints after the fix is confirmed.
5. If the root cause is genuinely unclear after investigation, **ask the user** rather than silencing the error.

---

## K562 Reference Gene Set

The following 10 protein-coding genes are the **preferred reference set** for simulations, analyses, figures, and PRO-seq processing across all tools (POLYARIS, KINEMA, GENEFIT, SCAN, ProtoPos). They are all expressed in K562 cells (human CML erythroleukemia), span 61-99 kb (optimised for simulation speed), and have 4-31 introns. GENCODE hg38 (GRCh38) canonical transcripts. Each gene has **measured elongation rates** from three independent DRB/Triptolide wavefront experiments (K562 BruDRB-seq, HEK293 ChIP-seq, HEK293 GRO-seq), enabling constrained parameter fitting in GENEFIT.

| # | Gene | Ensembl ID | Chr | Strand | Length (kb) | Introns | K562 Rate (kb/min) | CV | K562 Expr (RPKM) | Function |
|---|------|-----------|-----|--------|------------|---------|-------------------|-----|-----------------|----------|
| 1 | **USP10** | ENSG00000103194 | 16 | + | 80 | 13 | 1.80 | 0.031 | 2.61 | Deubiquitinase (p53 stabilisation). |
| 2 | **KIF2A** | ENSG00000068796 | 5 | + | 85 | 20 | 1.70 | 0.098 | 3.64 | Kinesin-13 motor protein (microtubule depolymerisation). |
| 3 | **RAP1B** | ENSG00000127314 | 12 | + | 62 | 7 | 1.84 | 0.113 | 4.06 | Ras-related small GTPase (integrin signalling). |
| 4 | **HNRNPC** | ENSG00000092199 | 14 | - | 62 | 8 | 0.98 | 0.124 | 6.66 | hnRNP C RNA-binding protein (mRNA processing). |
| 5 | **BAZ1B** | ENSG00000009954 | 7 | - | 82 | 19 | 1.62 | 0.137 | 2.77 | Bromodomain adjacent to zinc finger (chromatin remodeling). |
| 6 | **PDZD8** | ENSG00000165650 | 10 | - | 99 | 4 | 2.11 | 0.165 | 2.16 | PDZ domain-containing ER-mitochondria tether. |
| 7 | **ACSL3** | ENSG00000123983 | 2 | + | 84 | 16 | 2.49 | 0.192 | 3.39 | Long-chain acyl-CoA synthetase (lipid metabolism). |
| 8 | **CLIC4** | ENSG00000169504 | 1 | + | 99 | 5 | 2.23 | 0.203 | 2.37 | Chloride intracellular channel (redox signalling). |
| 9 | **CLTC** | ENSG00000141367 | 17 | + | 77 | 31 | 1.58 | 0.212 | 2.69 | Clathrin heavy chain (endocytosis). |
| 10 | **PPP2R2A** | ENSG00000221914 | 8 | + | 81 | 9 | 2.06 | 0.213 | 2.53 | PP2A regulatory subunit (cell cycle phosphatase). |

**When to use this gene set:**
- POLYARIS: Gene data files to be generated in `parameter_data/gene_data_files/`. Production configs to be updated.
- SCAN: Use these genes as the default analysis targets for figures and statistical comparisons.
- ProtoPos: Use these genes as the default targets for PRO-seq signal extraction and alignment QC.
- GENEFIT: Use these genes as fitting targets when optimising parameters against K562 data. Fix elongation rate from measured K562 BruDRB-seq values.

**Elongation rate sources:**
- K562 BruDRB-seq: Veloso et al. (2014) GSE55534 (SCAN experiment 2).
- HEK293 ChIP-seq: Cortazar et al. (2019) GSE134198 (SCAN experiment 4).
- HEK293 GRO-seq: Sheridan et al. (2019) GSE120201 (SCAN experiment 5).

**PPP termination probability (q) sources:** Genome-wide calibrated from Zimmer et al. (2021) PMID:34520723, Steurer et al. (2018) PMID:29632207, Gressel et al. (2019) PMID:31399571, Akcan & Heinig (2023) PMID:36727445. See `SCAN/data/k562_q_values.tsv`.

**Full qualifying gene catalogue:** 51 genes passing all quality filters are catalogued in SCAN experiment 6 (`SCAN/experiments/experiment6_reference_gene_selection/tables/all_qualifying_genes.tsv`). The 10 reference genes were selected from this catalogue as the lowest-CV gene per chromosome.

**Selection criteria:** 61-99 kb gene body (optimised for POLYARIS simulation speed), 4-31 introns (varied), K562 expression >= 2 RPKM, cross-validated elongation rates from 3 independent experiments (CV < 0.22), 10 different chromosomes. Covers deubiquitination, cytoskeleton, signalling, RNA processing, chromatin, ER-mitochondria, lipid metabolism, ion channels, endocytosis, and phosphatase regulation.

---

## Multi-Instance Coordination Protocol (v3)

Multiple Claude Code + Opus instances often run concurrently across the five repositories. This section defines the tiered permission model, master_agent role, quality gates, and autonomous operation rules.

**Full specification:** `trasim-lab/agents/multi_agent_work_plan.md` (v3).
**Living Spec:** `trasim-lab/configs/living_spec.md`.
**Interface contracts:** `trasim-lab/configs/interface_schema.yaml`.
**Status board:** `trasim-lab/agents/status_board.md`.

### Four-Tier Permission Model

| Tier | Name | Scope | Handoff? | Review? |
|------|------|-------|----------|---------|
| **0** | Read | Read any file in any repo. No writes outside own repo. | No | No |
| **1** | Soft Sync | Write to own repo + `trasim-lab/agents/`. Update shared docs. | No | No |
| **2** | Coordinated Integration | Cross-repo caller/callee changes. Both sides must update. | Yes (declared) | Daily reviewer |
| **3** | Structural / Contract Change | Modify `trasim-lab/configs/interface_schema.yaml`, rename/delete cross-repo symbols. | Yes + schema bump | Escalation reviewer |

### Per-Agent Write Access

| Instance launched from | Tier 0 | Tier 1 writes | Tier 2+ |
|------------------------|--------|---------------|---------|
| `POLYARIS/` | All repos | `POLYARIS/`, `trasim-lab/agents/` | Via handoff |
| `KINEMA/` | All repos | `KINEMA/`, `trasim-lab/agents/` | Via handoff |
| `GENEFIT/` | All repos | `GENEFIT/`, `trasim-lab/agents/` | Via handoff |
| `SCAN/` | All repos | `SCAN/`, `trasim-lab/agents/` | Via handoff |
| `ProtoPos/` | All repos | `ProtoPos/`, `trasim-lab/agents/` | Via handoff |
| **master_agent** (project root or trasim-lab) | All repos | All repos | All repos (declared) |

**Shared file ownership:**
- `GLOSSARY.md`: Owned by SCAN instance. Others read only; request SCAN to propagate updates via handoff.
- `trasim-lab/configs/interface_schema.yaml`: Any instance may propose changes but must bump the version and file handoffs to all consumers (Tier 3).
- `trasim-lab/configs/living_spec.md`: Human-owned. Agents may propose updates via handoff to SCAN.
- `trasim-lab/MASTER_CLAUDE.md`, `trasim-lab/CLAUDE.md`, root `environment.yml`, `setup_environment.sh`: Human-operated only.

### The master_agent

A Claude instance launched from the **project root** or **trasim-lab/** with cross-repo write access. Used for integration work, contract updates, and multi-repo refactors. Review rules: Tier 0-1 no review, Tier 2 daily review, Tier 3 escalation review. Must file an Integration Impact Declaration for any Tier 2+ session.

### Task IDs and Lifecycle

Every unit of work gets a Task ID: `T-<REPO>-<YYYYMMDD>-<NN>`. Handoffs get: `H-<SOURCE>-<TARGET>-<YYYYMMDD>-<NN>`. The master_agent uses `T-MASTER-<YYYYMMDD>-<NN>`.

**Task lifecycle:** `OPEN` -> `IN_PROGRESS` -> `REVIEW` -> `DONE`.
**Handoff lifecycle:** `NEW` -> `ACKNOWLEDGED` -> `IN_PROGRESS` -> `DONE`.

### Definition of Done

A task is not DONE unless **all** of these are true:
1. Code complete.
2. `pytest tests/ -v` passes with 0 failures. New logic has new tests.
3. At least one concrete example output (printed metrics, generated file, or screenshot).
4. Evidence bundle in the checkpoint (E1 diff, E2 repro commands, E3 test results, E4 example output, E5 assumptions).
5. Handoff filed if the change affects another repo.
6. Reviewer approved (`APPROVE` or `APPROVE WITH NITS`).

### Autonomous Operation Mode

When running unattended (user has launched the instance and walked away):

- **Do NOT use plan mode.** It requires user approval and blocks execution. Instead, write a brief plan to the checkpoint report and proceed.
- **Question budget:** Up to 3 clarification questions in the first 5 minutes of a session. After that, zero interactive questions. Make the conservative choice and document the decision in the checkpoint.
- **Never ask permission-style questions** ("Shall I proceed?", "Is this okay?"). Proceed with the best approach and document it.
- If genuinely uncertain about a destructive or irreversible action, **skip that subtask** and note it in the checkpoint rather than blocking the session.

### Checkpoint and Reporting Protocol

Every instance **must** write checkpoint reports to `trasim-lab/agents/checkpoints/`:

- **File name:** `<REPO>_checkpoint_<YYYYMMDD>_<HHMMSS>.md`
- **When:** Session start, every ~50 tool calls, session end, and on significant errors.
- **Required contents:** Task ID, status, **tier level**, files modified, evidence bundle (all 5 items), hallucination self-check, severity-tagged errors (`BLOCKER`/`WARNING`/`SUSPECT`), handoffs required.
- **Update `trasim-lab/agents/status_board.md`** at every checkpoint write.

### Two-Pass Internal Reviewer Loop

Reviews are internal to the Opus plan. No external model dependency.

- **Daily reviewer (Haiku):** Runs on every checkpoint at `REVIEW` status. Checks evidence completeness, hallucination self-check, obvious bugs, reproducibility, missing handoffs. For Tier 2: checks Integration Impact Declaration.
- **Escalation reviewer (Opus):** Triggered for Tier 3 changes, interface contracts, math/Markov models, results/figures, normalization, or when the daily reviewer flags issues. Checks scientific validity, units, coordinate conventions, result plausibility, Architecture Refactor Report.
- **Decisions:** `APPROVE` | `APPROVE WITH NITS` | `REQUEST CHANGES` | `INSUFFICIENT EVIDENCE`.

Run the collector: `bash trasim-lab/workflows/collect_reports.sh --daily` (or `--escalation`).

### Interface Contracts

`trasim-lab/configs/interface_schema.yaml` defines cross-repo function signatures, parameter types, units, and coordinate conventions. Versioned with MAJOR.MINOR.PATCH semver. Before modifying any provider function listed in a contract:
1. Check `trasim-lab/configs/interface_schema.yaml`.
2. If parameters/returns/types/units change: determine severity (PATCH/MINOR/MAJOR), bump `schema_version`, add `migration_log` entry, file handoffs to all consumers.

### Nightly Smoke Test

`bash trasim-lab/workflows/nightly_smoke_test.sh` runs basic imports, minimal tests, contract validation, Living Spec checks, and a toy cross-repo pipeline. Output: `trasim-lab/logs/smoke_test_<YYYYMMDD>.log` with per-repo PASS/FAIL.

### Git Commit Message Style

Commit messages must be **full descriptive sentences ending with a period**. They describe what was done, not imperative commands.

**Good examples:**
- `K562 10-gene reference dataset configured with gene-specific q values.`
- `Save simulation state added to the experiment.`
- `DRB and Triptolide fixed and now working from setup.`
- `Wire Cython-accelerated functions into TASEP hot loop for 2-3x speedup.`

**Bad examples:**
- `fix bug` (not a sentence, no period)
- `Update configs` (imperative, no period, not descriptive)
- `WIP` (not descriptive)

### Quality Control for Long Sessions

- Run `pytest tests/ -v` after every significant code change.
- If 3 or more consecutive test failures occur, write a checkpoint with severity `BLOCKER` and suggest the human start a fresh session.
- Never suppress errors with try/except workarounds.
- For sessions exceeding 4 hours of active coding, write a mid-session summary checkpoint.
- Complete the hallucination self-check honestly at every checkpoint.

---

## Research Direction Awareness

Claude agents have access to the project's research goals and execution plan in `trasim-lab/evaluation/research_goals.md` and `trasim-lab/evaluation/research_plan.md`. Agents should use this awareness to provide **opportunistic suggestions** about prompts, priorities, or approaches the user may not have considered.

### How It Works

This is a **zero-delay, end-of-session** mechanism. Agents must never pause their primary task to generate suggestions. The process is:

1. **During normal work**, the agent may notice things relevant to the research direction — a gap between what was asked and what the research goals call for, a missed opportunity, a dependency that should be addressed first, or a more efficient way to phrase a request. The agent does **not** stop to act on these observations.
2. **At the end of the session** (during the final checkpoint or after the last substantive task), the agent appends a brief entry to `trasim-lab/agents/research_suggestions.md`.
3. The user reviews this file asynchronously between sessions.

### When to Write a Suggestion

Write a suggestion **only** when one or more of these apply:

- The user's prompt targeted a subtask, but a different subtask would have unblocked more downstream work according to `trasim-lab/evaluation/research_plan.md`.
- The agent discovered during its work that a prerequisite for the user's goal is missing or incomplete (e.g., a dataset not yet acquired, a module not yet implemented).
- The user's prompt could have been more effective with a small change (e.g., specifying a gene name, a data type, or a particular phase from the research plan).
- The agent noticed an inconsistency between the current codebase state and the research goals (e.g., a claimed capability that doesn't exist yet, or a completed capability that the goals still list as pending).
- A cross-repo dependency or integration opportunity is ripe but hasn't been requested.

**Do not** write a suggestion if the session was straightforward and the user's prompt was well-targeted. Not every session needs a suggestion.

### Format

Append to `trasim-lab/agents/research_suggestions.md` using this format:

```markdown
### <YYYY-MM-DD> | <REPO> | <Task ID or "ad-hoc">

**Context:** <1 sentence: what the user asked for.>

**Suggestion:** <1-3 bullet points. Each bullet is a concrete, actionable suggestion. Reference specific research goals (e.g., "Goal 2"), phases (e.g., "Phase B"), genes, data types, or files where relevant.>

---
```

### Rules

- **Maximum 3 bullets per session.** Brevity is mandatory.
- **Never delay primary work.** If the session ends abruptly (e.g., context limit), skip the suggestion rather than sacrificing task completion.
- **Be concrete.** "Consider also fitting ChIP-seq data for USP10 (Goal 2, Phase B)" is useful. "Think about more data types" is not.
- **Reference the research documents.** Ground every suggestion in a specific goal, phase, risk, or success criterion from `trasim-lab/evaluation/research_goals.md` or `trasim-lab/evaluation/research_plan.md`.
- **No duplicates.** Before appending, scan the existing file to avoid repeating a suggestion already made.
- **Suggestions are advisory only.** The user decides whether to act on them. Never implement a suggestion without being asked.

---

## Interesting Results Protocol

Computational experiments sometimes produce findings that are scientifically interesting, unexpected, or that change the project's understanding of the system. These must be captured immediately in a standardised format so they are not lost between sessions.

### Directory and File Naming

All interesting results live in `trasim-lab/evaluation/interesting_results/`. Each result is a single Markdown file:

```
trasim-lab/evaluation/interesting_results/IR-<REPO>-<YYYYMMDD>-<NN>.md
```

Where `<REPO>` is the repository that produced the result (GENEFIT, KINEMA, POLYARIS, SCAN, ProtoPos), `<YYYYMMDD>` is the date, and `<NN>` is a zero-padded sequence number for that repo and date.

### Required Fields

Every interesting result file must contain the following fields:

```markdown
# IR-<REPO>-<YYYYMMDD>-<NN>: <Short Title>

**ID:** IR-<REPO>-<YYYYMMDD>-<NN>
**Date:** YYYY-MM-DD
**Status:** UNDOCUMENTED | DOCUMENTED
**Agent:** Which repo instance produced this
**Tool versions:** KINEMA commit, GENEFIT commit, etc.

---

## Intent
What experiment was being run and why.

## Method
Brief description of the computational method.

## Key Finding
1-3 sentence summary of the result.

## Evidence
Paths to CSV files, plots, logs.

## Interpretation
What this means biologically / for the project.

## Implications for trasim-lab/evaluation/
Which research plan items are affected.
```

### Status Lifecycle

- **UNDOCUMENTED**: The result has been recorded but not yet incorporated into `trasim-lab/evaluation/research_plan.md` or `trasim-lab/evaluation/research_goals.md`.
- **DOCUMENTED**: The nightly agent (or a human) has read the result, updated the evaluation documents accordingly, and marked it as processed.

### When to Create an Interesting Result

An agent should create an interesting result when any of these apply:

- An experiment produces a finding that contradicts a prior assumption in the research plan.
- A parameter degeneracy, identifiability issue, or unexpected invariant is discovered.
- A computational method fails in an informative way (e.g., demonstrates that a particular approach cannot work).
- A result confirms a hypothesis with quantitative evidence for the first time.
- A cross-repo integration reveals an unexpected interaction between models.

Routine successful fits, passing tests, and expected behaviour do NOT warrant an interesting result.

### Nightly Agent Behavior for UNDOCUMENTED Results

When the nightly smoke test flags UNDOCUMENTED results:

1. Read each UNDOCUMENTED result file.
2. Update `trasim-lab/evaluation/research_plan.md` accordingly (add findings, update timelines, note confirmed/refuted hypotheses).
3. Set status to DOCUMENTED in the result file.
4. Write a brief summary of what changed to the smoke test log.
5. Notify the human at the end of the nightly report.

-->
