// U1 convergence table: one row per quantity with pooled ESS and split R-hat and each Chain's ESS; sortable by ESS
// or R-hat, filterable by group and name. Clicking a row opens the trace panel.

import { getConvergence } from "../api";
import { filterQuantities, groups, quantities, sortQuantities, type Quantity, type SortKey } from "../convergence";
import { formatEss, formatRhat, RHAT_LIMIT } from "../format";
import type { RenderView } from "../view";
import { virtualList } from "../virtual-list";

const ROW_HEIGHT = 28;

function el<K extends keyof HTMLElementTagNameMap>(tag: K, props: Record<string, unknown> = {}, ...children: (Node | string)[]) {
  const node: HTMLElementTagNameMap[K] = Object.assign(document.createElement(tag), props);
  node.append(...children);
  return node;
}

const render: RenderView = (root, ctx) => {
  const chains = ctx.analysis.chains.map((c) => c.chain);
  const columns = `minmax(10em, 1.2fr) minmax(8em, 1.5fr) 6em 6em ${chains.map(() => "8em").join(" ")}`;

  const groupSelect = el("select", { title: "Group" }, el("option", { value: "", textContent: "All groups" }));
  const nameInput = el("input", { type: "search", placeholder: "Filter by name", value: ctx.params.get("name") ?? "" });
  const count = el("span", { className: "muted" });
  const header = el("div", { className: "table-row table-header" });
  header.style.gridTemplateColumns = columns;
  const body = el("div", { className: "table-body" });
  root.append(
    el("div", { className: "toolbar" }, groupSelect, nameInput, count),
    el("div", { className: "table" }, header, body),
  );

  let all: Quantity[] = [];
  let sortKey = (ctx.params.get("sort") as SortKey | null) ?? "table";
  let descending = ctx.params.get("desc") === "1";

  function headerCell(label: string, key: SortKey | null, title: string) {
    const cell = el("button", { className: "cell sortable", title, type: "button" }, label);
    if (key === null) {
      cell.disabled = true;
      return cell;
    }
    if (sortKey === key) cell.append(descending ? " ▼" : " ▲");
    cell.addEventListener("click", () => {
      // First click: the worst first (lowest ESS, highest R-hat); a second click reverses; a third restores the
      // table order.
      const worstFirst = key === "rhat";
      if (sortKey !== key) [sortKey, descending] = [key, worstFirst];
      else if (descending === worstFirst) descending = !descending;
      else [sortKey, descending] = ["table", false];
      update();
    });
    return cell;
  }

  function drawHeader() {
    header.replaceChildren(
      headerCell("Group", null, ""),
      headerCell("Name", null, ""),
      headerCell("ESS", "ess", "Multi-chain effective sample size; click to sort"),
      headerCell("R-hat", "rhat", "Split R-hat over the Chains; click to sort"),
      ...chains.map((c) => headerCell(`ESS ${c}`, null, `${c}'s own ESS`)),
    );
  }

  function renderRow(q: Quantity) {
    const flagged = q.rhat !== null && q.rhat > RHAT_LIMIT;
    const rhat = el("span", { className: `cell num${flagged ? " flag" : ""}` }, formatRhat(q.rhat) + (flagged ? " !" : ""));
    if (flagged) rhat.title = `R-hat above ${RHAT_LIMIT}: not converged`;
    const row = el(
      "div",
      { className: "table-row clickable", title: `Open the trace of ${q.group} ${q.name}` },
      el("span", { className: "cell muted" }, q.group),
      el("span", { className: "cell" }, q.name),
      el("span", { className: "cell num" }, formatEss(q.ess)),
      rhat,
      ...chains.map((c) => el("span", { className: "cell num muted" }, formatEss(q.chains[c] ?? null))),
    );
    row.style.gridTemplateColumns = columns;
    row.addEventListener("click", () => ctx.navigate("trace", { group: q.group, name: q.name }));
    return row;
  }

  const list = virtualList<Quantity>(body, ROW_HEIGHT, renderRow);

  function update() {
    const filter = { group: groupSelect.value, name: nameInput.value };
    const shown = sortQuantities(filterQuantities(all, filter), { key: sortKey, descending });
    count.textContent = `${shown.length.toLocaleString()} of ${all.length.toLocaleString()} quantities`;
    drawHeader();
    list.setItems(shown);
    const params: Record<string, string> = {};
    if (filter.group) params.group = filter.group;
    if (filter.name) params.name = filter.name;
    if (sortKey !== "table") Object.assign(params, { sort: sortKey, desc: descending ? "1" : "0" });
    ctx.replaceParams(params);
  }

  groupSelect.addEventListener("change", update);
  nameInput.addEventListener("input", update);
  count.textContent = "Loading the convergence table…";

  getConvergence()
    .then((rows) => {
      all = quantities(rows);
      groupSelect.append(...groups(all).map((g) => el("option", { value: g, textContent: g })));
      groupSelect.value = ctx.params.get("group") ?? "";
      update();
    })
    .catch((e: Error) => (count.textContent = `Could not load the convergence table: ${e.message}`));

  return () => list.destroy();
};

export default render;
