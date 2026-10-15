// U6 error table: an edge set's errors against the truth of a simulated input at a threshold. Missed true edges rank
// by how low their score is, false edges by how high; a false edge whose reverse is true is a direction error. Rows
// filter by kind and open the pair inspector, which shows the true Latent genotypes beside the posterior.

import { getJson } from "../api";
import { el, table } from "../dom";
import { formatRoute } from "../router";
import {
  ERROR_KINDS,
  errorCounts,
  errorsTable,
  filterErrors,
  matchText,
  pairLink,
  SOURCE_NAMES,
  type ErrorKind,
  type ErrorsPayload,
} from "../truth";
import type { RenderView } from "../view";

const fixed = (x: number | null, digits: number) => (x === null || Number.isNaN(x) ? "–" : x.toFixed(digits));
const KIND_LABEL = Object.fromEntries(ERROR_KINDS.map((k) => [k.kind, k.label])) as Record<ErrorKind, string>;

const render: RenderView = (root, ctx) => {
  let source = ctx.params.get("source") ?? "posterior";
  let threshold = Number(ctx.params.get("threshold") ?? 0.5);
  if (!(threshold >= 0 && threshold <= 1)) threshold = 0.5;
  const kindsParam = ctx.params.get("kinds");
  const kinds = new Set<ErrorKind>(
    kindsParam === null ? ERROR_KINDS.map((k) => k.kind) : (kindsParam.split(",").filter((k) => k in KIND_LABEL) as ErrorKind[]),
  );

  root.innerHTML = `
    <div class="toolbar">
      <h2>Errors against truth</h2>
      <label>Edges <select class="source"><option value="posterior">Posterior</option></select></label>
      <label>Called at score ≥ <input type="number" class="threshold" min="0" max="1" step="0.05" style="width: 5em"></label>
      <span class="kinds"></span>
    </div>
    <p class="truth muted"></p>
    <p class="status muted">Loading…</p>
    <div class="errors"></div>`;
  const $ = <T extends Element>(selector: string) => root.querySelector<T>(selector)!;
  const sourceSelect = $<HTMLSelectElement>(".source");
  const thresholdInput = $<HTMLInputElement>(".threshold");
  thresholdInput.value = String(threshold);
  let payload: ErrorsPayload | null = null;
  let request = 0;

  function saveParams() {
    ctx.replaceParams({ source, threshold: String(threshold), kinds: [...kinds].join(",") });
  }

  function drawKinds() {
    const counts = payload ? errorCounts(payload.rows) : null;
    $(".kinds").replaceChildren(
      ...ERROR_KINDS.map(({ kind, label, title }) => {
        const box = el("input", { type: "checkbox", checked: kinds.has(kind) });
        box.addEventListener("change", () => {
          if (box.checked) kinds.add(kind);
          else kinds.delete(kind);
          saveParams();
          drawTable();
        });
        return el("label", { title }, box, ` ${label}${counts ? ` (${counts[kind].toLocaleString()})` : ""}`);
      }),
    );
  }

  function drawTable() {
    if (!payload) return;
    const rows = filterErrors(payload.rows, kinds);
    const arrow = payload.directed ? "→" : "–";
    const directed = payload.directed;
    $(".status").textContent =
      `${rows.length.toLocaleString()} of ${payload.rows.length.toLocaleString()} errors of ${SOURCE_NAMES[payload.source]} ` +
      `(${payload.label}) at ${payload.threshold}, most confidently wrong first. Click a row to open the pair inspector.`;
    const head = ["Kind", "Edge", "Score", "Wrong by", ...(directed ? ["Reverse score", "Reverse true"] : []),
      "Lag", "True infection, parent", "True infection, child"];
    const numeric = head.flatMap((h, i) => (["Kind", "Edge", "Reverse true"].includes(h) ? [] : [i]));
    const t = table(head, rows.map((r) => [
      { content: KIND_LABEL[r.kind], className: r.kind === "direction" ? "flag" : "" },
      el("a", { href: pairLink(r) }, `${r.parent} ${arrow} ${r.child}`),
      fixed(r.marginal, 3),
      fixed(r.confidence, 3),
      ...(directed ? [fixed(r.reverse_marginal, 3), r.reverse_true ? "yes" : "no"] : []),
      fixed(r.lag, 1),
      fixed(r.parent_true_infection_time, 1),
      fixed(r.child_true_infection_time, 1),
    ]), { numeric });
    t.tBodies[0]!.querySelectorAll("tr").forEach((tr, i) => {
      for (const td of tr.children) td.classList.add("clickable");
      const r = rows[i]!;
      tr.addEventListener("click", (event) => {
        if ((event.target as Element).closest("a")) return;
        ctx.navigate("pair", { a: r.parent, b: r.child });
      });
    });
    $(".errors").replaceChildren(t);
  }

  function load() {
    const mine = ++request;
    $(".status").textContent = "Loading…";
    getJson<ErrorsPayload>("/api/errors", { source, threshold: String(threshold) })
      .then((p) => {
        if (mine !== request) return;
        payload = p;
        sourceSelect.replaceChildren(...p.sources.map((s) => el("option", { value: s, textContent: SOURCE_NAMES[s] })));
        sourceSelect.value = p.source;
        $(".truth").textContent = `${matchText(p.truth)}.`;
        drawKinds();
        drawTable();
      })
      .catch((e: Error) => {
        if (mine !== request) return;
        payload = null;
        $(".errors").replaceChildren();
        const status = $(".status");
        status.textContent = `No error table: ${e.message}. `;
        if (source !== "posterior") status.append(el("a", { href: formatRoute("errors") }, "Show the posterior's"));
      });
  }

  sourceSelect.addEventListener("change", () => {
    source = sourceSelect.value;
    saveParams();
    load();
  });
  thresholdInput.addEventListener("change", () => {
    const x = Number(thresholdInput.value);
    if (!(x >= 0 && x <= 1)) {
      thresholdInput.value = String(threshold);
      return;
    }
    threshold = x;
    saveParams();
    load();
  });

  drawKinds();
  load();

  // CSV export: the rows of the kinds shown, unrounded, under the payload's field names.
  return { table: () => (payload ? errorsTable(payload, filterErrors(payload.rows, kinds)) : null) };
};

export default render;
