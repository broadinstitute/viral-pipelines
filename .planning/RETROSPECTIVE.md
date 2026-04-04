# Project Retrospective

*A living document updated after each milestone. Lessons feed forward into future planning.*

---

## Milestone: v1.0 — VirNucPro Classification WDL

**Shipped:** 2026-04-01
**Phases:** 1 | **Plans:** 4 | **Tasks:** ~8

### What Was Built
- `classify_virnucpro_contigs` WDL task — inline Python3 heredoc with confidence-weighted delta scoring; classifies contigs as Viral/Non-viral/Ambiguous
- `classify_reads_by_contig` WDL task — inline DuckDB SQL pipeline joining PAF alignments to contig classifications; produces per-read call TSV
- Two standalone wrapper workflows with `allowNestedInputs: true` meta block, registered in `.dockstore.yml`
- Test input JSON files for both workflows; all WDL files pass `miniwdl check`

### What Worked
- **Selective `git show + sed` extraction:** Extracting only the new task lines (1488-1964) from the worktree branch avoided merging 30+ unrelated upstream changes while cleanly landing the new content on `ca-kb_python`
- **Inline Python heredocs:** The repo-wide `python3<<CODE` pattern passed miniwdl validation without any special handling; embedding scripts inline eliminated external file dependencies
- **Gap-closure plan:** Treating the worktree duplication as a discrete gap-closure plan (01-04) rather than a hotfix kept the audit trail clean and the verification cycle honest
- **DuckDB runtime pip install:** `pip install duckdb --quiet --no-cache-dir` before the Python block worked cleanly without a custom Docker image

### What Was Inefficient
- **Worktree content duplication:** Initial executor agents wrote content twice (lines duplicated verbatim), requiring a dedicated gap-closure plan. Root cause: agents appending to files without checking existing content.
- **Branch coordination overhead:** Plans 01–03 executed on separate worktree branches and required cherry-pick/merge steps before each plan. A single shared working branch would have eliminated this.
- **No `requirements-completed` in SUMMARY frontmatter for plans 01 and 03:** Made 3-source cross-reference rely entirely on VERIFICATION.md; plan 02 and 04 had the field, plans 01 and 03 did not.

### Patterns Established
- **WDL call alias pattern:** When a WDL workflow and task share the same name, use `call metagenomics.taskname as alias` on the call site. Applies to all future standalone wrapper workflows.
- **`allowNestedInputs: true` + workflow-level JSON keys:** Test input JSONs should use `{workflow_name}.{input_name}` format, not call-alias-qualified names.
- **Dockstore entries without `testParameterFiles` for placeholder JSONs:** Only add `testParameterFiles` when real runnable test data exists.

### Key Lessons
1. **Check file content before appending in executor agents.** Use `wc -l` or a content grep to verify the target file doesn't already contain the content before appending. This would have caught the duplication in plans 01–03.
2. **Single working branch preferred for multi-plan phases that all touch the same file.** When plans 01-03 all modify `tasks_metagenomics.wdl`, executing them on a shared branch avoids the cherry-pick/merge dance between plans.
3. **Gap-closure plans are the right pattern.** Rather than amending prior work or special-casing the executor, a clean gap-closure plan with explicit extraction steps was transparent, auditable, and easy to verify.

### Cost Observations
- Model mix: ~100% sonnet (all executor and verifier agents)
- Sessions: 1 (single session from plan to ship)
- Notable: Gap-closure plan added ~30% overhead to execution time but produced a cleaner audit trail than an in-place fix would have

---

## Milestone: v1.1 — Kraken2 Read Taxonomy Annotation WDL Task

**Shipped:** 2026-04-01
**Phases:** 1 | **Plans:** 2 | **Tasks:** 3

### What Was Built
- `parse_kraken2_reads` WDL task — inline DuckDB Python3 heredoc that joins Kraken2 per-read output to a pre-built taxonomy database; outputs 6-column TSV (SAMPLE_ID, READ_ID, TAXONOMY_ID, TAX_NAME, KINGDOM, TAX_RANK); 8 GB / 1 CPU / dynamic disk
- Standalone workflow `parse_kraken2_reads.wdl` using `as parse_reads` call alias
- Test input JSON with workflow-scoped placeholder keys; `.dockstore.yml` entry without `testParameterFiles`

### What Worked
- **DuckDB-only extraction discipline:** CONTEXT.md locked exactly which classes/functions to include vs exclude from the source script. The executor had unambiguous instructions and produced a clean heredoc on first attempt — no duplication, no rework.
- **Research phase confidence:** Researcher read the source Python script in full, catalogued the include/exclude list, and flagged the WDL Boolean-to-Python interpolation pitfall (`~{true="True" false="False" resolve_strains}`). Planner had high-confidence concrete values in every `<action>` block.
- **Worktree merge awareness:** Executor hit the old-base worktree issue (Wave 1 worktree lacked Phase 1 changes) and handled it automatically via `git merge ca-kb_python`. No manual intervention needed.
- **Single-task Phase 2:** Two plans with 1 task each made execution fast and verification tight (5/5 must-haves, zero gaps).

### What Was Inefficient
- **Repeated stash/merge conflicts on STATE.md:** Each worktree merge triggered a STATE.md conflict because the orchestrator's working tree had unstaged planning changes. Stashing before merge and popping after created repeated conflict cycles. Root cause: planning files modified in-session weren't committed before executor waves ran.
- **ROADMAP.md update lost to merge:** The `phase complete` CLI updated ROADMAP.md, but the subsequent stash-pop conflict resolution overwrote it. Required a manual re-edit at phase close.
- **`commit_docs: false` + git add for .planning:** The gitignore pattern `.planning/` (set by GSD for new files) meant some planning artifacts couldn't be staged with `git add .planning/`. Tracked files (ROADMAP.md, STATE.md) still worked; new files (VERIFICATION.md) were invisible to git. By design — but occasionally confusing.

### Patterns Established
- **Lock runtime parameters in CONTEXT.md:** When runtime specs are non-obvious (not matching the neighbor task), explicitly document them as locked decisions with values. Prevents executor guessing and plan checker warnings.
- **DuckDB-only extraction:** For tasks derived from multi-backend Python scripts, the CONTEXT.md `<decisions>` block should enumerate exactly which classes/functions to include and which to exclude. Concrete negatives ("do NOT include TaxonomyDatabase") are as important as positives.
- **Commit planning artifacts before executor waves:** To avoid STATE.md conflicts during worktree merges, stage and commit any in-session planning changes before spawning executor agents.

### Key Lessons
1. **Commit in-session planning changes before spawning executor worktrees.** Any uncommitted working-tree changes will create stash/merge conflicts when the worktree branch is merged back. A quick `git add -f .planning/ && git commit` before each wave prevents the conflict cycle.
2. **Phase complete CLI changes can be lost to merge.** After running `gsd-tools phase complete`, verify ROADMAP.md and STATE.md were actually updated on disk (not just in the CLI's memory). If a merge follows, re-check.
3. **High research confidence → zero executor rework.** v1.1 produced no duplication, no gap-closure plan, and a clean first-pass verification. The investment in the research phase (full script read, include/exclude list, pitfall documentation) paid off directly in execution quality.

### Cost Observations
- Model mix: ~100% sonnet (researcher, planner, checker, executor, verifier agents)
- Sessions: 1 (single session from plan-phase to complete-milestone)
- Notable: Zero gap-closure plans required (vs 1 in v1.0) — research phase investment paid off

---

## Cross-Milestone Trends

### Process Evolution

| Milestone | Phases | Plans | Key Change |
|-----------|--------|-------|------------|
| v1.0 | 1 | 4 | First milestone — baseline established |
| v1.1 | 1 | 2 | Tighter research phase → zero gap-closure plans |

### Cumulative Quality

| Milestone | miniwdl check | Verification | Gap-Closure Plans |
|-----------|--------------|--------------|-------------------|
| v1.0 | 3/3 pass | 11/11 (after gap closure) | 1 |
| v1.1 | 2/2 pass | 5/5 first pass | 0 |

### Top Lessons (Verified Across Milestones)

1. Content deduplication guard: always verify before appending to existing files
2. WDL standalone workflow call alias: required when task name equals workflow name
3. Commit in-session planning changes before executor worktree waves to avoid STATE.md merge conflicts
