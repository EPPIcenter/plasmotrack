// U6 error table: an edge set's errors against the truth of a simulated input at a threshold. Missed true edges rank
// by how low their score is, false edges by how high; a false edge whose reverse is true is a direction error. Rows
// filter by kind and open the pair inspector, which shows the true Latent genotypes beside the posterior.

import { getJson } from "../api";
import { formatRoute } from "../router";
import {
  ERROR_KINDS,
  errorCounts,
  filterErrors,
  matchText,
  pairLink,
  SOURCE_NAMES,
  type ErrorKind,
  type ErrorRow,
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
        const box = Object.assign(document.createElement("input"), { type: "checkbox", checked: kinds.has(kind) });
        box.addEventListener("change", () => {
          if (box.checked) kinds.add(kind);
          else kinds.delete(kind);
          saveParams();
          drawTable();
        });
        const l = Object.assign(document.createElement("label"), { title });
        l.append(box, ` ${label}${counts ? ` (${counts[kind].toLocaleString()})` : ""}`);
        return l;
      }),
    );
  }

  function cell(tag: "td" | "th", text: string, className = "") {
    return Object.assign(document.createElement(tag), { textContent: text, className });
  }

  function drawTable() {
    if (!payload) return;
    const rows = filterErrors(payload.rows, kinds);
    const arrow = payload.directed ? "→" : "–";
    $(".status").textContent =
      `${rows.length.toLocaleString()} of ${payload.rows.length.toLocaleString()} errors of ${SOURCE_NAMES[payload.source]} ` +
      `(${payload.label}) at ${payload.threshold}, most confidently wrong first. Click a row to open the pair inspector.`;
    const table = document.createElement("table");
    table.className = "data";
    const head = document.createElement("tr");
    const columns: [string, string][] = [
      ["Kind", ""], ["Edge", ""], ["Score", "num"], ["Wrong by", "num"],
      ...(payload.directed ? ([["Reverse score", "num"], ["Reverse true", ""]] as [string, string][]) : []),
      ["Lag", "num"], ["True infection, parent", "num"], ["True infection, child", "num"],
    ];
    head.append(...columns.map(([c, k]) => cell("th", c, k)));
    table.append(head);
    for (const r of rows) table.append(errorRow(r, arrow, payload.directed));
    $(".errors").replaceChildren(table);
  }

  function errorRow(r: ErrorRow, arrow: string, directed: boolean) {
    const tr = document.createElement("tr");
    const link = Object.assign(document.createElement("a"), { href: pairLink(r), textContent: `${r.parent} ${arrow} ${r.child}` });
    const edge = cell("td", "");
    edge.append(link);
    tr.append(
      cell("td", KIND_LABEL[r.kind], r.kind === "direction" ? "flag" : ""),
      edge,
      cell("td", fixed(r.marginal, 3), "num"),
      cell("td", fixed(r.confidence, 3), "num"),
      ...(directed ? [cell("td", fixed(r.reverse_marginal, 3), "num"), cell("td", r.reverse_true ? "yes" : "no")] : []),
      cell("td", fixed(r.lag, 1), "num"),
      cell("td", fixed(r.parent_true_infection_time, 1), "num"),
      cell("td", fixed(r.child_true_infection_time, 1), "num"),
    );
    for (const td of tr.children) td.classList.add("clickable");
    tr.addEventListener("click", (event) => {
      if ((event.target as Element).closest("a")) return;
      ctx.navigate("pair", { a: r.parent, b: r.child });
    });
    return tr;
  }

  function load() {
    const mine = ++request;
    $(".status").textContent = "Loading…";
    getJson<ErrorsPayload>("/api/errors", { source, threshold: String(threshold) })
      .then((p) => {
        if (mine !== request) return;
        payload = p;
        sourceSelect.replaceChildren(...p.sources.map((s) =>
          Object.assign(document.createElement("option"), { value: s, textContent: SOURCE_NAMES[s] })));
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
        if (source !== "posterior") status.append(Object.assign(document.createElement("a"), {
          href: formatRoute("errors"), textContent: "Show the posterior's",
        }));
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
};

export default render;
