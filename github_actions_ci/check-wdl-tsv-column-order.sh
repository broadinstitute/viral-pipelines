#!/bin/bash
# Guard the invariant that broke in PR #662: the assembly-metadata TSV column
# order must come from an explicit literal, never from the WDL keys() builtin.
#
# Cromwell evaluates a WDL map literal as `kvps.map(...).toMap`, producing a
# scala.collection.immutable.Map. That preserves insertion order only for 1-4
# entries; at 5+ it is a hash trie whose iteration order depends on the whole
# key set. So keys()/as_pairs() on the ~40-key stats_by_taxon maps return an
# arbitrary order on Terra, while miniwdl (whose Map is an ordered list of
# pairs) returns insertion order -- which is why no run-based test on either
# engine can catch this and the check has to be static.
#
# Symptom when this regresses: the entity:<table>_id column stops being column
# 1, Terra reads some other column as the entity id, and upload_entities_tsv
# fails with an opaque HTTP 400 "Duplicated entities are not allowed in TSV"
# -- or worse, succeeds with every column silently mis-populated.
#
#   ./github_actions_ci/check-wdl-tsv-column-order.sh

set -o pipefail

function absolute_path() {
    local SOURCE="$1"
    while [ -h "$SOURCE" ]; do
        DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
        if [[ "$OSTYPE" == "darwin"* ]]; then
            SOURCE="$(readlink "$SOURCE")"
        else
            SOURCE="$(readlink -f "$SOURCE")"
        fi
        [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
    done
    echo "$SOURCE"
}
SOURCE="${BASH_SOURCE[0]}"
SCRIPT=$(absolute_path "$SOURCE")
SCRIPT_DIRNAME="$(dirname "$SOURCE")"
SCRIPTPATH="$(cd -P "$(echo "$SCRIPT_DIRNAME")" &> /dev/null && pwd)"
REPO_PATH="$(realpath "$SCRIPTPATH/../")"

cd "${REPO_PATH}"

# This check parses WDL with miniwdl's own AST, so it needs an interpreter that
# can import WDL. In CI that is plain python3 (tests install miniwdl with pip3).
# Locally miniwdl is often a pipx/venv install invisible to the system python3,
# so fall back to the interpreter named in the miniwdl launcher's shebang.
PYTHON=""
for candidate in python3 python; do
    if command -v "$candidate" &> /dev/null && "$candidate" -c 'import WDL' &> /dev/null; then
        PYTHON="$candidate"
        break
    fi
done
if [ -z "$PYTHON" ] && command -v miniwdl &> /dev/null; then
    shebang_py="$(head -1 "$(command -v miniwdl)" | sed 's|^#!||' | awk '{print $1}')"
    if [ -x "$shebang_py" ] && "$shebang_py" -c 'import WDL' &> /dev/null; then
        PYTHON="$shebang_py"
    fi
fi
if [ -z "$PYTHON" ]; then
    echo "ERROR: could not find a python interpreter that can 'import WDL'." >&2
    echo "       install miniwdl (pip3 install miniwdl) and re-run." >&2
    exit 1
fi

echo "Checking WDL TSV column-order invariants (using $PYTHON)"

"$PYTHON" - <<'PYTHON_EOF'
import sys
import glob
import WDL

# (workflow, array-literal decl, map-literal decl) triples whose column order
# must agree. Both names must exist in the file; a rename must break this check
# loudly rather than silently disable it.
TARGETS = [
    ("pipes/WDL/workflows/assemble_denovo_metagenomic.wdl",   "assembly_header", "stats_by_taxon"),
    ("pipes/WDL/workflows/scaffold_and_refine_multitaxa.wdl", "assembly_header", "stats_by_taxon"),
]

# Ordering of these depends on the engine's Map iteration order, which Cromwell
# does not make insertion-ordered. as_map() is deliberately NOT banned: it is
# used legitimately for write_json() output, where key order is irrelevant.
BANNED_FUNCTIONS = ("keys", "as_pairs")

IMPORT_PATH = ["pipes/WDL/tasks"]
problems = []


def src_text(lines, pos):
    """Exact source text of an expression (miniwdl end_column is one past the end)."""
    if pos.line == pos.end_line:
        return lines[pos.line - 1][pos.column - 1 : pos.end_column - 1]
    out = [lines[pos.line - 1][pos.column - 1 :]]
    out += lines[pos.line : pos.end_line - 1]
    out.append(lines[pos.end_line - 1][: pos.end_column - 1])
    return "".join(out)


def column_name(lines, expr):
    """A constant string yields its value; one with a placeholder (e.g.
    "~{table_name}") has no literal, so compare its source text instead. Both
    sides of the comparison are spelled the same way, so this matches cleanly
    without having to evaluate anything."""
    literal = getattr(expr, "literal", None)
    if literal is not None:
        return literal.value
    return src_text(lines, expr.pos)


def walk(node):
    yield node
    for child in getattr(node, "children", []):
        yield from walk(child)


def local_nodes(doc):
    """Nodes defined in this document only -- not pulled in via imports, which
    would otherwise be reported once per importer."""
    if doc.workflow:
        yield from walk(doc.workflow)
    for task in doc.tasks:
        yield from walk(task)


# ---- 1. column order must match between the literal and the map ------------
for path, header_name, map_name in TARGETS:
    try:
        doc = WDL.load(path, path=IMPORT_PATH)
    except Exception as exc:
        problems.append("{}: failed to parse: {}".format(path, exc))
        continue
    lines = open(path).read().splitlines()

    header_decl = map_decl = None
    for node in local_nodes(doc):
        if isinstance(node, WDL.Tree.Decl):
            if node.name == header_name and isinstance(node.expr, WDL.Expr.Array):
                header_decl = node
            elif node.name == map_name and isinstance(node.expr, WDL.Expr.Map):
                map_decl = node

    if header_decl is None:
        problems.append(
            "{}: no Array literal declaration named '{}' found. If it was renamed or "
            "removed, update TARGETS in this script -- do not leave the check "
            "silently matching nothing.".format(path, header_name))
    if map_decl is None:
        problems.append(
            "{}: no Map literal declaration named '{}' found. If it was renamed or "
            "removed, update TARGETS in this script -- do not leave the check "
            "silently matching nothing.".format(path, map_name))
    if header_decl is None or map_decl is None:
        continue

    header_cols = [column_name(lines, e) for e in header_decl.expr.items]
    map_cols = [column_name(lines, k) for k, _ in map_decl.expr.items]

    # (a) Set equality is mandatory: the workflow scatters over the header array
    # and looks each name up in the map, and a missing Map key is a hard error.
    missing_from_map = [c for c in header_cols if c not in set(map_cols)]
    missing_from_header = [c for c in map_cols if c not in set(header_cols)]
    if missing_from_map:
        problems.append(
            "{}: {} (Ln {}) names column(s) absent from {} (Ln {}): {}. "
            "At runtime this is a 'Map key not found' failure.".format(
                path, header_name, header_decl.pos.line, map_name,
                map_decl.pos.line, missing_from_map))
    if missing_from_header:
        problems.append(
            "{}: {} (Ln {}) has key(s) absent from {} (Ln {}): {}. Those columns "
            "would be silently dropped from the emitted TSV.".format(
                path, map_name, map_decl.pos.line, header_name,
                header_decl.pos.line, missing_from_header))

    # (b) Order equality is not required for correctness once the header drives
    # the row extraction, but keeping the two lists in the same order is free and
    # keeps them reviewable side by side.
    if not missing_from_map and not missing_from_header and header_cols != map_cols:
        first = next(i for i, (a, b) in enumerate(zip(header_cols, map_cols)) if a != b)
        problems.append(
            "{}: {} (Ln {}) and {} (Ln {}) list the same columns in different "
            "orders; first difference at index {}: {!r} vs {!r}.".format(
                path, header_name, header_decl.pos.line, map_name,
                map_decl.pos.line, first, header_cols[first], map_cols[first]))

    # (c) No duplicate column names in the header.
    dups = sorted({c for c in header_cols if header_cols.count(c) > 1})
    if dups:
        problems.append(
            "{}: {} (Ln {}) repeats column name(s) {}.".format(
                path, header_name, header_decl.pos.line, dups))

    # (d) The Terra invariant that actually broke: Terra reads the entity id out
    # of column 1 regardless of that column's name.
    if not header_cols[0].startswith("entity:"):
        problems.append(
            "{}: {} (Ln {}) must start with the Terra entity id column "
            "('entity:<table>_id'), but column 1 is {!r}. Terra keys each row on "
            "column 1 whatever it is named.".format(
                path, header_name, header_decl.pos.line, header_cols[0]))

    print("  {}: {} n={}, {} n={}, order matches, entity id first".format(
        path, header_name, len(header_cols), map_name, len(map_cols)))

# ---- 2. no WDL-native keys()/as_pairs() anywhere in pipes/WDL --------------
# Matching on the AST rather than by grep avoids false positives on the many
# Python dict .keys() calls inside task command blocks.
for path in sorted(glob.glob("pipes/WDL/workflows/*.wdl") + glob.glob("pipes/WDL/tasks/*.wdl")):
    try:
        doc = WDL.load(path, path=IMPORT_PATH)
    except Exception as exc:
        problems.append("{}: failed to parse: {}".format(path, exc))
        continue
    for node in local_nodes(doc):
        if isinstance(node, WDL.Expr.Apply) and node.function_name in BANNED_FUNCTIONS:
            problems.append(
                "{} (Ln {}, Col {}): {}() is banned -- its result order depends on "
                "the engine's Map iteration order, which Cromwell does not make "
                "insertion-ordered. Use an explicit Array[String] literal for "
                "column order instead. See PR #662 and the notes atop this "
                "script.".format(path, node.pos.line, node.pos.column,
                                 node.function_name))

if problems:
    print("\nFAILED -- {} problem(s):".format(len(problems)), file=sys.stderr)
    for p in problems:
        print("  * {}".format(p), file=sys.stderr)
    sys.exit(1)

print("\nAll WDL TSV column-order invariants OK.")
PYTHON_EOF
