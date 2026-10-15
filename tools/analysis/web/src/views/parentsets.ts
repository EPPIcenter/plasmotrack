// U6 parent-set panel: per Infection, the effective parent-set size, the credible collection's size at a level,
// whether the true set is in it, and the modal set's probability, beside the null-model fit when there is one.
// Sortable by any column; a summary line per fit.

import { getJson } from "../api";
import { sortRows, summaryLine, type ParentSetKey, type ParentSetRow, type ParentSetsPayload } from "../parentsets";
import type { RenderView } from "../view";
import { virtualList } from "../virtual-list";

const ROW_HEIGHT = 28;

interface Column {
  key: ParentSetKey;
  label: string;
  title: string;
  format: (row: ParentSetRow) => string;
}

const number = (x: number | undefined, digits: number) => (x === undefined ? "–" : x.toFixed(digits));
const yesNo = (x: boolean | null | undefined) => (x === null || x === undefined ? "–" : x ? "yes" : "no");

function columns(withNull: boolean, withTruth: boolean): Column[] {
  const fit = (prefix: "" | "null_"): Column[] => {
    const label = (text: string) => (prefix ? `Null: ${text.toLowerCase()}` : text);
    return [
      { key: `${prefix}effective_size`, label: label("Effective size"), title: "exp of the entropy of the pooled Parent-set weights",
        format: (r) => number(r[`${prefix}effective_size`], 2) },
      { key: `${prefix}credible_size`, label: label("Credible size"), title: "Sets in the credible collection at the level",
        format: (r) => String(r[`${prefix}credible_size`] ?? "–") },
      ...(withTruth
        ? [{ key: `${prefix}covered` as ParentSetKey, label: label("True set in it"), title: "Whether the true Parent set is in the credible collection",
             format: (r: ParentSetRow) => yesNo(r[`${prefix}covered`]) }]
        : []),
      { key: `${prefix}top_probability`, label: label("Top probability"), title: "Pooled weight of the modal Parent set",
        format: (r) => number(r[`${prefix}top_probability`], 3) },
    ];
  };
  return [
    { key: "infection", label: "Infection", title: "Infection id", format: (r) => r.infection },
    ...fit(""),
    ...(withNull ? fit("null_") : []),
  ];
}

const render: RenderView = (root, ctx) => {
  const level = Number(ctx.params.get("level") ?? "0.8");
  let sortKey = (ctx.params.get("sort") as ParentSetKey | null) ?? "infection";
  let descending = ctx.params.get("desc") === "1";

  root.innerHTML = `
    <div class="toolbar">
      <h2>Parent sets</h2>
      <label>Credible level <input type="number" class="level" min="0.05" max="1" step="0.05"></label>
      <span class="count muted"></span>
    </div>
    <ul class="summary"></ul>
    <div class="table"><div class="table-row table-header"></div><div class="table-body"></div></div>`;
  const $ = <T extends Element>(selector: string) => root.querySelector<T>(selector)!;
  const levelInput = $<HTMLInputElement>(".level");
  levelInput.value = String(level);
  levelInput.addEventListener("change", () => {
    const value = Number(levelInput.value);
    if (value > 0 && value <= 1) ctx.navigate("parentsets", { level: String(value) });
  });
  const header = $<HTMLElement>(".table-header");
  $(".count").textContent = "Loading…";

  let payload: ParentSetsPayload | null = null;
  let cols: Column[] = [];
  let template = "";

  function renderRow(r: ParentSetRow) {
    const row = document.createElement("div");
    row.className = "table-row ruled";
    row.style.gridTemplateColumns = template;
    for (const c of cols) {
      const missed = (c.key === "covered" || c.key === "null_covered") && r[c.key] === false;
      row.append(Object.assign(document.createElement("span"), {
        className: `cell${c.key === "infection" ? "" : " num"}${missed ? " flag" : ""}`,
        textContent: c.format(r),
      }));
    }
    return row;
  }
  const list = virtualList<ParentSetRow>($(".table-body"), ROW_HEIGHT, renderRow);

  function update() {
    if (!payload) return;
    header.style.gridTemplateColumns = template;
    header.replaceChildren(
      ...cols.map((c) => {
        const cell = Object.assign(document.createElement("button"), {
          className: "cell sortable", type: "button", title: `${c.title}; click to sort`, textContent: c.label,
        });
        if (sortKey === c.key) cell.append(descending ? " ▼" : " ▲");
        cell.addEventListener("click", () => {
          if (sortKey === c.key) descending = !descending;
          else [sortKey, descending] = [c.key, c.key !== "infection"];
          update();
        });
        return cell;
      }),
    );
    list.setItems(sortRows(payload.rows, sortKey, descending));
    ctx.replaceParams({ level: String(level), ...(sortKey !== "infection" || descending ? { sort: sortKey, desc: descending ? "1" : "0" } : {}) });
  }

  getJson<ParentSetsPayload>("/api/parent-sets", { level: String(level) })
    .then((p) => {
      payload = p;
      const withNull = p.fits.some((f) => f.fit === "null");
      cols = columns(withNull, p.truth);
      template = cols.map((c) => (c.key === "infection" ? "minmax(6em, 1fr)" : "minmax(7em, 1fr)")).join(" ");
      if (!cols.some((c) => c.key === sortKey)) sortKey = "infection";
      $(".count").textContent = `${p.rows.length.toLocaleString()} Infections` + (p.truth ? "" : "; no truth, so no coverage");
      $(".summary").replaceChildren(...p.fits.map((f) => Object.assign(document.createElement("li"), { textContent: summaryLine(f, p.level) })));
      update();
    })
    .catch((e: Error) => ($(".count").textContent = `Could not load the parent sets: ${e.message}`));

  return () => list.destroy();
};

export default render;
