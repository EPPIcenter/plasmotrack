// U10 infection card: what the posterior says about one Infection. Densities per Chain of the expected false positive
// and false negative counts, duration and infection time; COI per Draw; the Latent genotype as a Locus by Allele
// heatmap of posterior presence with the Observed Alleles outlined; the Parent-set distribution with the credible
// collection marked; parent and child Edge marginals linking to the pair inspector and the trace panel; q, entropy and
// effective size. Linked as `#/infection?id=<id>`.

import * as Plot from "@observablehq/plot";
import { getJson, num } from "../api";
import { CHAIN_COLOURS, formatEss } from "../format";
import {
  coiBars, density, heatmapCells, infectionParams, readInfection, setLabel,
  type EdgeMarginal, type InfectionResponse, type Quantity, type SetRow,
} from "../infection";
import { pairParams } from "../pair";
import { formatRoute } from "../router";
import type { RenderView } from "../view";

const INK = "#52514e";
const GRID = "#d6d5d0";
// Posterior presence: the reference palette's sequential blue, from near the surface (absent) to its darkest step.
const PRESENCE = { type: "linear" as const, domain: [0, 0.5, 1], range: ["#f1f0ec", "#6da7ec", "#0d366b"], label: "Posterior presence" };
const CREDIBLE = "#2a78d6";
const OUTSIDE = "#b8b7b1";

function el<K extends keyof HTMLElementTagNameMap>(tag: K, props: Record<string, unknown> = {}, ...children: (Node | string)[]) {
  const node: HTMLElementTagNameMap[K] = Object.assign(document.createElement(tag), props);
  node.append(...children);
  return node;
}

/** A probability to three decimals, or to two significant digits below 0.01 so small weights are not rounded away. */
const prob = (x: number | null) => (x === null ? "–" : x > 0 && x < 0.01 ? x.toPrecision(2) : x.toFixed(3));
const fixed = (x: number | null, digits = 3) => (x === null ? "–" : x.toFixed(digits));

function table(head: string[], rows: (string | Node)[][]) {
  return el(
    "table",
    { className: "data" },
    el("thead", {}, el("tr", {}, ...head.map((h) => el("th", {}, h)))),
    el("tbody", {}, ...rows.map((cells) => el("tr", {}, ...cells.map((c) => el("td", {}, c))))),
  );
}

const QUANTITIES: { key: keyof InfectionResponse["quantities"]; title: string; digits: number }[] = [
  { key: "false_positives", title: "Expected false positives", digits: 4 },
  { key: "false_negatives", title: "Expected false negatives", digits: 4 },
  { key: "duration", title: "Duration", digits: 1 },
  { key: "infection_time", title: "Infection time", digits: 1 },
];

const render: RenderView = (root, ctx) => {
  const id = readInfection(ctx.params);
  const chains = ctx.analysis.chains.map((c) => c.chain);
  const colour = { domain: chains, range: CHAIN_COLOURS.slice(0, chains.length) };
  const card = (other: string) => el("a", { href: formatRoute("infection", infectionParams(other)) }, other);

  const input = el("input", { type: "search", placeholder: "Infection id", value: id, size: 12 });
  const form = el("form", { className: "toolbar" }, el("h2", {}, "Infection"), input, el("button", { type: "submit" }, "Open"));
  const body = el("div", { className: "infection" });
  root.append(form, body);
  form.addEventListener("submit", (e) => {
    e.preventDefault();
    ctx.navigate("infection", infectionParams(input.value));
  });

  if (!id) {
    body.append(el("p", { className: "muted" }, "Enter an Infection id to see its posterior."));
    return;
  }
  body.textContent = "Loading…";

  getJson<InfectionResponse>("/api/infection", { id })
    .then((p) => {
      body.replaceChildren();
      const section = (title: string, ...children: (Node | string)[]) => body.append(el("h3", {}, title), ...children);
      const s = p.summary;

      // Summary: q, expected observed parents, entropy, effective size
      section(`Infection ${p.infection}, observed at ${fixed(p.observation_time, 1)}`,
        table(["Chain", "q", "Expected observed parents", "Entropy", "Effective size"], [
          ...s.chains.map((c) => [c.chain, prob(c.q), fixed(c.expected_observed_parents), fixed(c.entropy), fixed(c.effective_size)]),
          ["All Chains", prob(s.q), fixed(s.expected_observed_parents), fixed(s.entropy), fixed(s.effective_size)],
        ]),
        el("p", { className: "muted" }, "Means over Draws: q, the weight of Parent sets holding the unobserved source; ",
          "entropy of the Parent-set weights and its exponential, the effective number of Parent sets. ",
          el("a", { href: formatRoute("trace", { group: "unobserved-source inclusion", name: p.infection }) }, "Trace of q"), "."));

      // Densities
      const densities = el("div", { className: "histograms" });
      for (const { key, title, digits } of QUANTITIES) {
        const q: Quantity = p.quantities[key];
        const curves = q.chains.flatMap((c) => density(c.value.map((v) => num(v) as number), 120).map((d) => ({ chain: c.chain, ...d })));
        const caption = `median ${fixed(q.median, digits)}, 95% ${fixed(q.lower, digits)} – ${fixed(q.upper, digits)}`;
        const heading = q.group ? el("a", { href: formatRoute("trace", { group: q.group, name: q.name ?? "" }) }, title) : title;
        densities.append(el("figure", {}, el("figcaption", {}, heading, el("span", { className: "muted" }, ` · ${caption}`)),
          Plot.plot({
            width: 320,
            height: 160,
            color: colour,
            x: { label: null, grid: true },
            y: { label: "Density", ticks: 3 },
            marks: [
              Plot.ruleX([q.lower, q.upper].filter((v) => v !== null), { stroke: INK, strokeDasharray: "3,3" }),
              Plot.ruleX([q.median].filter((v) => v !== null), { stroke: INK }),
              Plot.lineY(curves, { x: "x", y: "y", z: "chain", stroke: "chain", strokeWidth: 2 }),
              Plot.tip(curves, Plot.pointerX({ x: "x", y: "y", channels: { Chain: "chain" } })),
              Plot.ruleY([0], { stroke: GRID }),
            ],
          })));
      }
      section("Observation process, duration and infection time",
        el("div", { className: "legend" }, Plot.legend({ color: { ...colour, type: "categorical" } })),
        densities,
        el("p", { className: "muted" }, "Kernel densities per Chain; solid rule: pooled median; dashed: pooled 2.5% and 97.5% ",
          "quantiles. Infection time is the observation time minus the duration."));

      // COI
      const coiRows = coiBars(p.coi.chains);
      section("COI",
        Plot.plot({
          width: 320 + 40 * p.coi.distribution.length,
          height: 180,
          color: colour,
          fx: { label: "COI (largest Allele count at any Locus of the Latent genotype)" },
          x: { axis: null, domain: chains },
          y: { label: "Draws", grid: true },
          marks: [
            Plot.barY(coiRows, { fx: "coi", x: "chain", y: "draws", fill: "chain", inset: 1, tip: true }),
            Plot.ruleY([0], { stroke: GRID }),
          ],
        }),
        el("p", { className: "muted" }, "All Chains: ", p.coi.distribution.map((d) => `COI ${d.coi} in ${(100 * d.fraction).toFixed(1)}% of Draws`).join("; "), "."));

      // Genotype heatmap
      const cells = heatmapCells(p.genotype.loci);
      const maxAlleles = Math.max(0, ...p.genotype.loci.map((l) => l.alleles));
      section("Latent genotype",
        Plot.plot({
          width: Math.max(320, 70 + 28 * maxAlleles),
          height: 40 + 28 * p.genotype.loci.length,
          marginLeft: 70,
          padding: 0,
          color: { ...PRESENCE, legend: true },
          x: { label: "Allele", domain: Array.from({ length: maxAlleles }, (_, k) => k) },
          y: { label: null, domain: p.genotype.loci.map((l) => l.locus) },
          marks: [
            Plot.cell(cells, { x: "allele", y: "locus", fill: "presence", inset: 3, tip: true,
              channels: { Observed: (c: (typeof cells)[number]) => (c.observed ? "yes" : "no") } }),
            Plot.cell(cells.filter((c) => c.observed), { x: "allele", y: "locus", fill: "none", stroke: "#0b0b0b", strokeWidth: 2, inset: 1 }),
          ],
        }),
        el("p", { className: "muted" }, "Fill: fraction of Draws of every Chain whose Latent genotype has the Allele. Ringed: ",
          "Alleles in the Observed genotype; a Locus with none outlined is unobserved."));

      // Parent sets
      const level = p.parent_sets.level;
      const setPlot = (title: string, rows: SetRow[]) => {
        const shown = rows.map((r) => ({ ...r, label: `${setLabel(r.members, p.infection)} ${r.class}` }));
        return el("figure", {}, el("figcaption", {}, title, el("span", { className: "muted" },
          ` · credible collection: ${rows.filter((r) => r.credible).length} of ${rows.length}`)),
        Plot.plot({
          width: 520,
          height: 30 + 24 * shown.length,
          marginLeft: 260,
          x: { label: "Pooled weight", domain: [0, 1], grid: true },
          y: { label: null, domain: shown.map((r) => r.label) },
          color: { domain: ["in the credible collection", "outside it"], range: [CREDIBLE, OUTSIDE] },
          marks: [
            Plot.barX(shown, { x: "weight", y: "label", insetTop: 2, insetBottom: 2,
              fill: (r: SetRow) => (r.credible ? "in the credible collection" : "outside it"),
              tip: true, channels: { Rank: "rank", "Cumulative mass": "cumulative_mass", "Set id": "set_id" } }),
            Plot.text(shown, { x: "weight", y: "label", text: (r: SetRow) => prob(r.weight), textAnchor: "start", dx: 4, fill: INK }),
            Plot.ruleX([0], { stroke: GRID }),
          ],
        }));
      };
      section("Parent sets",
        el("div", { className: "legend" }, Plot.legend({ color: { domain: ["in the credible collection", "outside it"], range: [CREDIBLE, OUTSIDE] } })),
        el("div", { className: "histograms" }, setPlot("All Chains", p.parent_sets.pooled),
          ...p.parent_sets.chains.map((c) => setPlot(c.chain, c.sets))),
        el("p", { className: "muted" }, `Sets by pooled weight (mean over Draws), with composition class. Credible collection: `,
          `the smallest prefix with mass at least ${level}. All Chains weights each Chain by its Draws and matches sets by members.`));

      // Edge marginals
      const edgeList = (title: string, rows: EdgeMarginal[], other: (e: EdgeMarginal) => string) => {
        if (!rows.length) return el("p", { className: "muted" }, "None.");
        return table([title, "Pooled", "", ...chains, "ESS", ""], rows.map((e) => {
          const who = other(e);
          const bar = el("span", { className: "bar" });
          bar.style.width = `${Math.round(100 * e.pooled)}%`;
          const links = el("span", {}, el("a", { href: formatRoute("trace", { group: e.group, name: e.name }) }, "trace"));
          if (!e.unobserved) links.append(" · ", el("a", { href: formatRoute("pair", pairParams({ a: e.parent, b: e.child })) }, "pair"));
          return [e.unobserved ? el("span", { className: "muted" }, "unobserved source") : card(who),
            el("span", { className: "bar-track" }, bar), prob(e.pooled),
            ...chains.map((ch) => prob(e.chains.find((c) => c.chain === ch)?.mean ?? 0)), formatEss(num(e.ess)), links];
        }));
      };
      section(`Parents of ${p.infection}`, edgeList("Parent", p.parents, (e) => e.parent));
      section(`Children of ${p.infection}`, edgeList("Child", p.children, (e) => e.child));
      body.append(el("p", { className: "muted" }, "Edge marginals: the mean over Draws of the summed weight of the child's Parent sets holding the parent."));
    })
    .catch((e: Error) => (body.textContent = `Could not load Infection ${id}: ${e.message}`));
};

export default render;
