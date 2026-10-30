---
name: miint-review
description: Run the miint-reviewer agent to convergence on the current branch — one /code-review pass, then review, fix, re-review — with a disposition ledger and a stop rule. Use whenever a review of pending changes is asked for: "miint-reviewer", "miint-style review", "review this branch", "review my changes", "review before I open the PR". Owns the loop only; the rules live in the agent.
---

# Reviewing a branch to convergence

The first review pass finds what reading the diff finds;
the fixes then change the diff, and the next pass reads what the fixes wrote. This
skill owns that loop. The rules, the trigger map, and the report format live in
`.claude/agents/miint-reviewer.md` — do not restate them here or in your prompts.

**You are the caller.** The reviewer is read-only by design: its instructions rule
out building, running tests or probes, and editing files. Everything it hands back
that requires *doing* — applying a fix, rebuilding, running a test file, running a
`[needs-probe]` — is yours.

The base is the branch the PR targets. It is named for the DuckDB release
(`v1.5-variegata`), not `main`.

## Before round 1

**`/simplify` is optional, and runs at most once, first.** It edits code rather than
reporting, so it cannot go through the ledger, and run between rounds it rewrites code
a reviewer has already read. Its reuse and simplification lenses overlap the
reviewer's. Run it when the branch adds a lot of new code, where duplication is
likely; skip it for a focused fix. Its edits are then part of the diff everything
after reads.

**`/code-review` runs once, before the reviewer rounds, and is not run again.** Invoke
the `code-review` skill at `medium` effort on the branch, report-only: no `--fix`, no
`--comment`. If the Skill tool refuses it, ask the user to type
`/code-review medium`. Every finding it reports goes into the ledger with a
disposition, like a reviewer finding. It runs first because bug fixes rewrite the
most code, and the reviewer rounds that follow then read code close to its final
shape and check the fixes as part of their ask. Re-running it on a changed diff
re-raises what the ledger already settled.

## Every round gets its own reviewer

**Spawn a fresh `miint-reviewer` for each round.** A round is a verdict on the code as
it stands now, and a reviewer that already concluded something about a span does not
re-derive it — it remembers. Remembering is the failure mode: it clears code it
cleared before the fix rewrote it, and it defends its earlier findings instead of
re-testing them against what is there. The check that a fix worked has to come from a
reader that did not watch it being made.

Independence is the default; you steer through the **prompt**, not through the agent's
memory. What you carry across rounds is the ledger, not the reviewer.

**`SendMessage` is for steering a reviewer within its own round**, where its context is
the point — it continues the same agent, which still holds the diff it read and the
spans it opened. Use it to:

- hand back probe results it asked for, so it re-tags `[needs-probe]` → `[verified]`
  or drops the finding;
- ask what a finding means, or which of two readings of a rule it applied;
- narrow it mid-flight ("skip R7 on the vendored file, it isn't ours").

Never use it to start a new round. If the code changed, that is a new reviewer.

## The loop

**Round 1.**

```
Agent(subagent_type: "miint-reviewer",
      prompt: "Review the current branch against <base>. <anything branch-specific>")
```

Include the `/code-review` findings and their dispositions in the prompt, so the
reviewer does not re-raise what you already settled.

**Between rounds — disposition, then fix.**

Give every entry under the report's **Findings** exactly one disposition, and write it
into the ledger:

| Disposition | Means | Next |
|---|---|---|
| `fixed` | You changed the code | Name the commit or `path:line` |
| `declined` | Correct rule, wrong here | Record the **reason** — it goes into the next prompt |
| `deferred` | Real, but wider than this PR (the reviewer's Scope discipline) | Record the issue's number, or its draft title and body until it is opened (see the PR-body block), then drop it from the loop |
| `probed` | Was `[needs-probe]`; you ran it | Record the result — confirmed or refuted |

**Questions** get an answer, recorded once in the ledger; a later round does not
re-raise an answered question. **Follow-ups** are pre-existing problems: they need no
disposition, and they go into the closing prose list, not the loop. **Unprobed** lines
go into the closing prose list too, and into the PR block's `Not probed` when they bear
on this PR.

Then apply the fixes and run what the change touches: rebuild (`bash build.sh`), each
changed or affected SQL test file (`./build/release/test/unittest test/sql/<file>.test`),
the C++ tests for a changed module (`./build/release/extension/miint/tests "[<tag>]"`, or
`tests_duckdb` for the tests that open a DuckDB),
and `make format-check` in the pinned environment (`CLAUDE.md`, "Formatting"). Run
`MIINT_TEST_JOBS=8 bash run_tests.sh` (its SQL and C++ tests on eight workers) before the last
round. Fix test failures before the next round, so it reads code that passes. Commit
the fixes before the next round: a reviewer reads
`git diff <base>...HEAD`, so a fix left in the working tree is invisible to it. A test
the handoff lists as skipped in CI is run here, or recorded under "Not probed" with the
reason.

**Round N+1 — a new agent, carrying the ledger in its prompt.** Three things, and
nothing else:

1. **Scope** — the paths that changed since the last round, and any rule family
   already settled. This is what keeps a fresh spawn from costing a full-branch read.
2. **The standing ledger** — every `declined` and `deferred` finding and every answered
   question, with its reason. This is what stops round 3 from re-raising what round 2
   settled.
3. **The ask** — review the current state of the named paths, confirm each `fixed`
   finding is actually fixed, and flag anything the fixes introduced.

State in the prompt that the ledger is the *caller's* disposition, not a finding the
reviewer made — it did not make these calls and should not treat them as its own.

## When to stop

There is no round cap. Stop when any of these holds, and say which one in your summary:

- **Nothing new.** Every finding in this round is one already in the ledger.
- **Nothing actionable.** The new findings are all ones you decline, on reasons already
  in the ledger — the reviewer is re-litigating a class you have settled.
- **Churn.** The round's findings are mostly about code the previous round's fixes
  wrote, and the fixes were right. Two passes disagreeing about a rewrite is a signal
  to stop and let a human read it, not to write a third version.
- **What's left isn't a fix.** Everything open is `deferred` (an issue) or
  `needs-probe` (run the probe, or file it — do not spend a round on it).
- **Diminishing returns.** The findings are real but minor, and each round is
  returning less than the last. Say so and stop.

A stop is not a pass. End with what remains open and why, so the human reviewer
inherits the ledger rather than re-deriving it.

## Cost discipline

The reviewer has its own, in the agent file. Yours:

- **Scope the prompt.** A fresh reviewer re-reads the whole branch unless the prompt
  names a slice. Round 2+ names the changed paths and the rules still in play, so it
  reads only those.
- **Carry the ledger.** A decline the reviewer can see is a decline it does not
  re-derive.
- **Keep each reviewer short-lived.** A reviewer carried across rounds accumulates every
  diff and span it has read. Per-round spawns keep each one's context to a single round.
- **Don't fan out from here.** The reviewer decides its own fan-out. A second layer of
  subagents on top duplicates its work.

## Report each round

Keep it to a table plus the stop decision. The user is watching convergence, not
re-reading findings:

```
code-review — 5 findings: 3 fixed, 1 declined (test-only path), 1 probed (refuted)
Round 1 — 9 findings, 4 questions: 6 fixed, 2 declined, 1 deferred (#N); questions answered
Round 2 — 2 findings: both repeats of declined. Stopping: nothing actionable.
Open: 1 deferred (#N), 0 needs-probe.
```

## Always end with the deferred list, in prose

**When you stop, the last thing you write is a plain-language list of everything left
open** — every `deferred` finding, every unrun `needs-probe`, every follow-up the
reviewer proposed, and anything you declined that a human might reasonably re-decide.
It is separate from the per-round table above: the table tracks convergence, this
hands over the work.

The reader is deciding what to do next, so write it for someone who has not read the
diff:

- **Group by what the reader would do about it** — cleanups that are somebody's next
  PR, probes that need an environment you lacked (a WASM build, macOS, the live
  service), and things only they can settle (a modelling call, a product decision,
  filing an issue you have no access to file).
- **One line each: what it is, where, and why it is still open.** Name the file, not
  the finding number — `R6-3` means nothing tomorrow. "Wider than this PR" and
  "needs a WASM build" are reasons; "deferred" is not.
- **Say what a probe would settle**, not just that it is unrun.
- **Separate what is pre-existing from what this branch introduced.** A reviewer
  inheriting the list needs to know which items the PR is answerable for.

Keep it short — a few grouped lines, not a report. If nothing is open, say that in one
sentence rather than omitting the section, so the reader can tell the difference
between a clean stop and a forgotten one.

## Then write the PR-body block

A PR records that this loop ran and what it left open. `review-loop-check` in CI fails
a PR whose description lacks the block (the `no-agent-review` label opts out; see
`CLAUDE.md`, "Reviewing a branch"). After the prose list, print the block filled in from
the ledger (or, while a deferred finding has no issue number, what the last rule below
says instead). It replaces the template's whole `## Reviewer loop` section, heading
included:

```markdown
## Reviewer loop

- Reviewed at: <short sha of the last commit a round read> · rounds: <N> · code-review: <level> · stopped: <which stop rule>
- Fixed: <count>
- Declined:
  - <what, `path`> — <why>
- Deferred:
  - <what, `path`> — #<issue>
- Not probed:
  - <the question> — <why it was not run>
```

- `Reviewed at` is a commit on the branch. If commits land after the last round, say so
  on that line ("plus 1 later commit: format fix") instead of moving the sha. A rebase
  rewrites the reviewed commit: name its rebased counterpart, and say so on the line
  ("rebased from abc1234 onto the new base").
- Write `none` under a heading with no entries; do not drop the heading.
- A deferred item names its issue by number. Every deferred issue is opened before
  the PR is, so the description carries `#N` from its first version; `issue to file`
  is not an entry. Open them yourself only when the author has asked you to. Otherwise,
  or when you lack access, put each issue's title and body after the prose list in
  place of the block, and print the block once the author gives you the numbers.
- Do not edit the PR description yourself unless asked; hand the block to the author.
