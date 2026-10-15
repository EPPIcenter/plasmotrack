// U3 pair inspector: for Infections A and B, P(A before B) per Chain with switch counts and the per-Draw order trace;
// both Edge marginals per Chain; co-inclusion of A and B in third Infections' Parent sets; genotypes side by side per
// Locus with differing Loci highlighted; posterior infection times and the lag. Linked as `#/pair?a=<id>&b=<id>`.

import * as Plot from "@observablehq/plot";
import { getJson } from "../api";
import { el, table } from "../dom";
import { CHAIN_COLOURS } from "../format";
import { differingLoci, pairParams, pairTables, readPair, type GenotypeRow, type PairResponse } from "../pair";
import { formatRoute } from "../router";
import type { RenderView } from "../view";

const INK = "#52514e";

/** A probability to three decimals, or to two significant digits below 0.01 so small weights are not rounded away. */
const prob = (x: number | null) => (x === null ? "–" : x > 0 && x < 0.01 ? x.toPrecision(2) : x.toFixed(3));
const time = (x: number | null) => (x === null ? "–" : x.toFixed(1));

function genotype(bits: string | null | undefined) {
  return bits ? el("code", {}, bits) : el("span", { className: "muted" }, "unobserved");
}

const render: RenderView = (root, ctx) => {
  const pair = readPair(ctx.params);
  const chains = ctx.analysis.chains.map((c) => c.chain);
  const colour = { domain: chains, range: CHAIN_COLOURS.slice(0, chains.length) };

  const inputA = el("input", { type: "search", placeholder: "Infection A", value: pair.a, size: 12 });
  const inputB = el("input", { type: "search", placeholder: "Infection B", value: pair.b, size: 12 });
  const swap = el("button", { type: "button", title: "Swap A and B" }, "⇄");
  const form = el("form", { className: "toolbar" }, el("h2", {}, "Pair inspector"), "A ", inputA, swap, "B ", inputB,
    el("button", { type: "submit" }, "Inspect"));
  const body = el("div", { className: "pair" });
  root.append(form, body);
  form.addEventListener("submit", (e) => {
    e.preventDefault();
    ctx.navigate("pair", pairParams({ a: inputA.value, b: inputB.value }));
  });
  swap.addEventListener("click", () => ctx.navigate("pair", pairParams({ a: inputB.value, b: inputA.value })));

  if (!pair.complete) {
    body.append(el("p", { className: "muted" }, pair.a && pair.a === pair.b
      ? "Enter two different Infection ids."
      : "Enter two Infection ids to see their order, edges, genotypes and infection times."));
    return;
  }
  const { a, b } = pair;
  body.textContent = "Loading…";
  let drawTrace = () => {};
  let loaded: PairResponse | null = null;

  getJson<PairResponse>("/api/pair", { a, b })
    .then((p) => {
      loaded = p;
      body.replaceChildren();
      const section = (title: string, ...children: (Node | string)[]) => body.append(el("h3", {}, title), ...children);

      // Order
      const traceBox = el("div", { className: "trace-plot" });
      section(`Order: P(${a} before ${b})`,
        table(["Chain", "Draws", `P(${a} before ${b})`, "Switches"], [
          ...p.order.chains.map((c) => [c.chain, String(c.draws), prob(c.probability), String(c.switches)]),
          ["All Chains", String(p.order.chains.reduce((n, c) => n + c.draws, 0)), prob(p.order.probability), ""],
        ]),
        el("p", { className: "muted" }, "Switches: adjacent Draws at which the order changes. The trace is 1 at Draws where ",
          `${a} precedes ${b} in the Ordering.`),
        el("div", { className: "legend" }, Plot.legend({ color: { ...colour, type: "categorical" } })),
        traceBox);
      drawTrace = () => {
        const points = p.order.chains.flatMap((c) => c.iteration.map((iteration, i) => ({ chain: c.chain, iteration, before: c.before[i]! })));
        traceBox.replaceChildren(Plot.plot({
          width: Math.max(480, traceBox.clientWidth),
          height: 70 * chains.length + 40,
          marginLeft: 70,
          color: colour,
          x: { label: "Sampling iteration", grid: true },
          y: { domain: [0, 1], inset: 8, ticks: [0, 1], tickFormat: (d: number) => (d ? `${a} first` : `${b} first`), label: null },
          fy: { label: null, domain: chains, padding: 0.3 },
          marks: [
            Plot.frame({ stroke: "#d6d5d0" }),
            Plot.lineY(points, { x: "iteration", y: "before", fy: "chain", stroke: "chain", strokeWidth: 2, curve: "step-after" }),
            Plot.tip(points, Plot.pointerX({ x: "iteration", y: "before", fy: "chain",
              channels: { Chain: "chain" }, format: { y: (d: number) => (d ? `${a} first` : `${b} first`), fy: false } })),
          ],
        }));
      };
      drawTrace();

      // Edge marginals
      const edgeLink = (child: string, parent: string, pooled: number | null) =>
        pooled ? el("a", { href: formatRoute("trace", { group: "edge", name: `${child}<-${parent}` }) }, `${child} ← ${parent}`)
          : `${child} ← ${parent}`;
      section("Edge marginals",
        table(["Chain", `${b} ← ${a}`, `${a} ← ${b}`], [
          ...p.edges.chains.map((c) => [c.chain, prob(c.a_to_b), prob(c.b_to_a)]),
          ["All Chains", prob(p.edges.a_to_b), prob(p.edges.b_to_a)],
        ]),
        el("p", { className: "muted" }, "Mean over Draws of the summed weight of the child's Parent sets holding the parent. Traces: ",
          edgeLink(b, a, p.edges.a_to_b), ", ", edgeLink(a, b, p.edges.b_to_a), "."));

      // Co-inclusion
      section(`Co-inclusion of ${a} and ${b} in a third Infection's Parent sets`,
        p.co_inclusion.length
          ? table(["Infection", "All Chains", ...chains], p.co_inclusion.map((r) =>
            [r.infection, prob(r.weight), ...chains.map((ch) => prob(r.chains.find((c) => c.chain === ch)?.weight ?? 0))]))
          : el("p", { className: "muted" }, `No Infection has a Parent set holding both ${a} and ${b}.`));

      // Genotypes
      const rows: GenotypeRow[] = p.genotypes.rows;
      const differing = differingLoci(rows);
      const truth = p.genotypes.truth;
      section("Genotypes",
        el("p", {}, el("span", { className: differing.count ? "flag" : "" }, `${differing.text}.`), " ",
          el("span", { className: "muted" }, "Latent: posterior mode over every Chain's Draws, with the fraction of Draws that had it.")),
        table([
          "Locus", `Observed ${a}`, `Latent ${a}`, "agreement", `Latent ${b}`, "agreement", `Observed ${b}`,
          ...(truth ? [`True ${a}`, `True ${b}`] : []),
        ], rows.map((r) => [
          r.locus, genotype(r.observed_a), genotype(r.latent_a), prob(r.agreement_a), genotype(r.latent_b), prob(r.agreement_b),
          genotype(r.observed_b), ...(truth ? [genotype(r.true_a), genotype(r.true_b)] : []),
        ]), { rowClass: (i) => (rows[i]!.differs ? "differs" : "") }));

      // Infection times
      const times = [{ infection: a, ...p.times.a }, { infection: b, ...p.times.b }];
      const intervals = times.filter((t) => t.lower !== null);
      section("Infection times",
        table(["Infection", "Observation time", "Infection time, median", "95% interval"], times.map((t) =>
          [el("a", { href: `#/infection?id=${encodeURIComponent(t.infection)}` }, t.infection), time(t.observation_time), time(t.median), t.lower === null ? "–" : `${time(t.lower)} – ${time(t.upper)}`])),
        el("p", {}, `Lag (observation time of ${b} minus ${a}): `, el("strong", {}, time(p.times.lag)), "."),
        Plot.plot({
          width: 640,
          height: 110,
          marginLeft: 70,
          x: { label: "Time", grid: true },
          y: { domain: [a, b], label: null },
          marks: [
            Plot.link(intervals, { x1: "lower", x2: "upper", y1: "infection", y2: "infection", stroke: INK, strokeWidth: 2 }),
            Plot.dot(intervals, { x: "median", y: "infection", fill: INK, r: 4, tip: true,
              title: (t: typeof times[number]) => `${t.infection}: median ${time(t.median)}, 95% ${time(t.lower)} – ${time(t.upper)}` }),
            Plot.dot(times, { x: "observation_time", y: "infection", symbol: "times", stroke: "#0b0b0b", r: 5, strokeWidth: 2 }),
          ],
        }),
        el("p", { className: "muted" }, "Line: 95% posterior interval of infection time (observation time minus duration quantiles); ",
          "dot: median; ×: observation time."));
    })
    .catch((e: Error) => (body.textContent = `Could not load the pair: ${e.message}`));

  const resize = new ResizeObserver(() => drawTrace());
  resize.observe(body);
  return {
    cleanup: () => resize.disconnect(),
    // CSV export: each section's numbers under the payload's field names, one file per section.
    tables: () => (loaded ? pairTables(loaded) : null),
  };
};

export default render;
