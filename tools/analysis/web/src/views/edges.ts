// U2 edge table: one row per (child, parent) with the pooled Edge marginal, between-Chain spread, switch count, ESS
// and each Chain's mean and switch count; sortable by every column, filterable by child or parent id, with the
// unobserved source shown or hidden. Clicking a row opens the edge's trace panel.

import { getJson } from "../api";
import { edgeRows, filterEdges, sortEdges, type EdgeJson, type EdgeRow, type SortKey } from "../edges";
import { formatEss } from "../format";
import { getSelection, subscribeSelection } from "../selection";
import type { RenderView } from "../view";
import { virtualList } from "../virtual-list";

const ROW_HEIGHT = 28;

function el<K extends keyof HTMLElementTagNameMap>(tag: K, props: Record<string, unknown> = {}, ...children: (Node | string)[]) {
  const node: HTMLElementTagNameMap[K] = Object.assign(document.createElement(tag), props);
  node.append(...children);
  return node;
}

const formatProb = (x: number | undefined) => (x === undefined ? "–" : x.toFixed(3));

/** The direction of a column's first click: the rows worth looking at first (disagreeing, probable, locked, poorly
 * mixed). */
function worstFirst(key: SortKey): boolean {
  return key === "spread" || key === "pooled" || key.startsWith("mean:");
}

const render: RenderView = (root, ctx) => {
  const chains = ctx.analysis.chains.map((c) => c.chain);
  const columns = `minmax(6em, 1fr) minmax(10em, 1.4fr) 6em 6em 6em 6em ${chains.map(() => "8.5em 8.5em").join(" ")}`;

  const childInput = el("input", { type: "search", placeholder: "Child id", value: ctx.params.get("child") ?? "" });
  const parentInput = el("input", { type: "search", placeholder: "Parent id", value: ctx.params.get("parent") ?? "" });
  const unobserved = el("input", { type: "checkbox", checked: ctx.params.get("unobserved") !== "0" });
  // The Infections selected in the network view (#84): when there are any, the table keeps their edges unless unticked.
  const selected = el("input", { type: "checkbox", checked: true });
  const selectedLabel = el("label", {}, selected, " Selected Infections only");
  const count = el("span", { className: "muted" });
  const header = el("div", { className: "table-row table-header" });
  header.style.gridTemplateColumns = columns;
  const body = el("div", { className: "table-body" });
  root.append(
    el("div", { className: "toolbar" }, childInput, parentInput, el("label", {}, unobserved, " Unobserved source"), selectedLabel, count),
    el("div", { className: "table" }, header, body),
  );

  let all: EdgeRow[] = [];
  let sortKey = (ctx.params.get("sort") as SortKey | null) ?? "table";
  let descending = ctx.params.get("desc") === "1";

  function headerCell(label: string, key: SortKey, title: string) {
    const cell = el("button", { className: "cell sortable", title: `${title}; click to sort`, type: "button" }, label);
    if (sortKey === key) cell.append(descending ? " ▼" : " ▲");
    cell.addEventListener("click", () => {
      // First click: the rows worth looking at first; a second click reverses; a third restores the table order.
      const first = worstFirst(key);
      if (sortKey !== key) [sortKey, descending] = [key, first];
      else if (descending === first) descending = !descending;
      else [sortKey, descending] = ["table", false];
      update();
    });
    return cell;
  }

  function drawHeader() {
    header.replaceChildren(
      headerCell("Child", "child", "The Infection whose parents these are"),
      headerCell("Parent", "parent", "An observed parent, or the unobserved source (the child's Latent parent)"),
      headerCell("Pooled", "pooled", "Edge marginal, the mean over every Draw of every Chain"),
      headerCell("Spread", "spread", "Largest minus smallest Chain mean; a Chain without the edge counts as 0"),
      headerCell("Switches", "switches", "Adjacent Draws at which the marginal crossed one half, summed over Chains"),
      headerCell("ESS", "ess", "Multi-chain effective sample size"),
      ...chains.flatMap((c) => [
        headerCell(`Mean ${c}`, `mean:${c}`, `${c}'s mean over its Draws`),
        headerCell(`Switches ${c}`, `switches:${c}`, `${c}'s crossings of one half`),
      ]),
    );
  }

  function renderRow(r: EdgeRow) {
    const parent = r.unobserved
      ? el("span", { className: "cell muted", title: `${r.parent}, the unobserved source` }, "unobserved source")
      : el("span", { className: "cell" }, r.parent);
    const row = el(
      "div",
      { className: "table-row clickable", title: `Open the trace of ${r.group} ${r.name}` },
      el("span", { className: "cell" }, r.child),
      parent,
      el("span", { className: "cell num" }, formatProb(r.pooled)),
      el("span", { className: "cell num" }, formatProb(r.spread)),
      el("span", { className: "cell num" }, String(r.switches)),
      el("span", { className: "cell num" }, formatEss(r.ess)),
      ...chains.flatMap((c) => [
        el("span", { className: "cell num muted" }, formatProb(r.chains[c]?.mean)),
        el("span", { className: "cell num muted" }, r.chains[c] ? String(r.chains[c].switches) : "–"),
      ]),
    );
    row.style.gridTemplateColumns = columns;
    row.addEventListener("click", () => ctx.navigate("trace", { group: r.group, name: r.name }));
    return row;
  }

  const list = virtualList<EdgeRow>(body, ROW_HEIGHT, renderRow);

  function update() {
    const selection = getSelection();
    selectedLabel.hidden = selection.size === 0;
    selectedLabel.lastChild!.textContent = ` Selected Infections only (${[...selection].slice(0, 5).join(", ")}${selection.size > 5 ? ", …" : ""})`;
    const filter = { child: childInput.value, parent: parentInput.value, unobserved: unobserved.checked,
                     selection: selected.checked ? selection : undefined };
    const shown = sortEdges(filterEdges(all, filter), { key: sortKey, descending });
    count.textContent = `${shown.length.toLocaleString()} of ${all.length.toLocaleString()} edges`;
    drawHeader();
    list.setItems(shown);
    const params: Record<string, string> = {};
    if (filter.child.trim()) params.child = filter.child.trim();
    if (filter.parent.trim()) params.parent = filter.parent.trim();
    if (!filter.unobserved) params.unobserved = "0";
    if (sortKey !== "table") Object.assign(params, { sort: sortKey, desc: descending ? "1" : "0" });
    ctx.replaceParams(params);
  }

  childInput.addEventListener("input", update);
  parentInput.addEventListener("input", update);
  unobserved.addEventListener("change", update);
  selected.addEventListener("change", update);
  const unsubscribe = subscribeSelection(() => all.length && update());
  count.textContent = "Loading the edge table…";

  getJson<{ rows: EdgeJson[] }>("/api/edges")
    .then(({ rows }) => {
      all = edgeRows(rows);
      update();
    })
    .catch((e: Error) => (count.textContent = `Could not load the edge table: ${e.message}`));

  return () => {
    unsubscribe();
    list.destroy();
  };
};

export default render;
