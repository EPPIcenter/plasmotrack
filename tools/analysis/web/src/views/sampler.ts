// U4 sampler panel, per Chain: acceptance rate per sampler over Windows (a small multiple per sampler type, a line per
// target, the target band shaded), swap acceptance per temperature pair, the temperature ladder over adaptation, the
// burn-in-end proposal variances, and the Run events as rules on every time axis.

import * as Plot from "@observablehq/plot";
import { getJson } from "../api";
import { el, table } from "../dom";
import { CHAIN_COLOURS } from "../format";
import {
  describeEvent,
  EVENT_NAMES,
  eventRules,
  groupByType,
  ladderPoints,
  parseBand,
  severityOf,
  timeDomain,
  varianceRows,
  windowX,
  type EventRule,
  type PanelWindow,
  type SamplerPanel,
  type Severity,
} from "../sampler";
import type { RenderView } from "../view";

const BURN_IN_GREY = "#a3a29d";
const MUTED = "#52514e";
const BAND = "#2a78d6";
/** Status colours of the reference palette; each rule also names its events in its tooltip and the legend above. */
const SEVERITY_COLOUR: Record<Severity, string> = { critical: "#d03b3b", warning: "#fab219", neutral: MUTED };
/** More series than categorical colours: thin grey lines, named in the tooltip, never cycled hues. */
const MAX_COLOURED = CHAIN_COLOURS.length;

const STYLE = `
  .sampler-panel .multiples { display: flex; flex-wrap: wrap; gap: 1.2em; }
  .sampler-panel .multiples figure { margin: 0; }
  .sampler-panel input[type="number"] { width: 4.5em; }
  .sampler-panel table { border-collapse: collapse; font-variant-numeric: tabular-nums; }
  .sampler-panel th, .sampler-panel td { text-align: left; padding: 0.2em 0.8em 0.2em 0; border-bottom: 1px solid var(--rule); vertical-align: top; }
  .sampler-panel td.num { text-align: right; }
  .sampler-panel .event-key { display: inline-flex; align-items: center; gap: 0.3em; margin-right: 1em; }
  .sampler-panel .event-key i { display: inline-block; width: 0; height: 1em; border-left: 2px solid; }
  .sampler-panel .event-key i.neutral { border-left-style: dashed; }
`;

const fmt = (x: number) => Number(x.toPrecision(3)).toString();
const range = (w: PanelWindow) => `Window ${w.window} (${w.phase.replace("_", "-")}, iterations ${w.first_iteration}–${w.last_iteration})`;

const render: RenderView = (root, ctx) => {
  const chains = ctx.analysis.chains.map((c) => c.chain);
  let chain = ctx.params.get("chain") ?? chains[0] ?? "";
  let band = parseBand(ctx.params.get("band_lo"), ctx.params.get("band_hi"));
  let panel: SamplerPanel | null = null;
  let request = 0;

  root.innerHTML = `
    <style>${STYLE}</style>
    <div class="sampler-panel">
      <div class="toolbar">
        <h2>Samplers</h2>
        <label>Chain <select class="chain"></select></label>
        <label>Target acceptance band <input type="number" class="band-lo" min="0" max="1" step="0.05">
          to <input type="number" class="band-hi" min="0" max="1" step="0.05"></label>
      </div>
      <p class="events-key"></p>
      <h3>Acceptance rate per Window <span class="muted">(accepted over proposed, Hot chain; a small multiple per
        sampler type, a line per target; shaded: the target band; grey: burn-in)</span></h3>
      <div class="multiples acceptance"></div>
      <h3>Swap acceptance per temperature pair</h3>
      <div class="swaps"></div>
      <h3>Temperature ladder <span class="muted">(at each adaptation during burn-in, then the final ladder)</span></h3>
      <div class="ladder"></div>
      <h3 class="variances-title">Final proposal variances</h3>
      <div class="variances"></div>
      <h3>Run events</h3>
      <div class="events"></div>
    </div>`;
  const $ = <T extends Element>(selector: string) => root.querySelector<T>(selector)!;

  const select = $<HTMLSelectElement>(".chain");
  select.append(...chains.map((c) => el("option", { value: c, textContent: c })));
  select.value = chain;
  const [loInput, hiInput] = [$<HTMLInputElement>(".band-lo"), $<HTMLInputElement>(".band-hi")];
  const showBand = () => ([loInput.value, hiInput.value] = [String(band[0]), String(band[1])]);
  showBand();

  const keepParams = () =>
    ctx.replaceParams({ chain, band_lo: String(band[0]), band_hi: String(band[1]) });
  select.addEventListener("change", () => {
    chain = select.value;
    keepParams();
    load();
  });
  for (const input of [loInput, hiInput]) {
    input.addEventListener("change", () => {
      band = parseBand(loInput.value, hiInput.value);
      showBand();
      keepParams();
      drawAcceptance();
    });
  }

  /** Burn-in shading and the event rules, shared by every time axis so the rules line up across plots. */
  function timeMarks(p: SamplerPanel, rules: EventRule[]): Plot.Markish[] {
    const burnIn = p.windows.filter((w) => w.phase === "burn_in");
    const marks: Plot.Markish[] = [];
    if (burnIn.length) {
      const [x1, x2] = timeDomain(burnIn);
      marks.push(Plot.rect([{ x1, x2 }], { x1: "x1", x2: "x2", fill: BURN_IN_GREY, fillOpacity: 0.12 }));
    }
    const rule = (dashed: boolean) =>
      Plot.ruleX(rules.filter((r) => (r.severity === "neutral") === dashed), {
        x: "iteration",
        stroke: (r: EventRule) => SEVERITY_COLOUR[r.severity],
        strokeWidth: 2,
        strokeDasharray: dashed ? "4 3" : undefined,
        title: "label",
      });
    marks.push(
      rule(true),
      rule(false),
      Plot.dot(rules, {
        x: "iteration",
        frameAnchor: "top",
        symbol: "triangle",
        rotate: 180,
        r: 5,
        fill: (r: EventRule) => SEVERITY_COLOUR[r.severity],
        title: "label",
      }),
    );
    return marks;
  }

  /** A colour scale over `series` when they fit the categorical palette; null for more, drawn grey. */
  function colours(series: string[]) {
    return series.length <= MAX_COLOURED
      ? { domain: series, range: CHAIN_COLOURS.slice(0, series.length), legend: series.length > 1 }
      : null;
  }

  function drawAcceptance() {
    if (!panel) return;
    const p = panel;
    const domain = timeDomain(p.windows);
    const rules = eventRules(p.events);
    const groups = groupByType(p.acceptance);
    if (!groups.length) {
      $(".acceptance").textContent = "No acceptance rows in the closed Windows.";
      return;
    }
    $(".acceptance").replaceChildren(
      ...groups.map((g) => {
        const data = g.rows.map((r) => ({ ...r, x: windowX(r) }));
        const colour = colours(g.targets);
        const stroke = colour ? "target" : MUTED;
        return Plot.plot({
          width: 380,
          height: 200,
          marginLeft: 44,
          marginTop: 24,
          title: g.type,
          subtitle:
            `${g.targets.length} target${g.targets.length === 1 ? "" : "s"}` +
            (g.rows.every((r) => r.rate === null) ? ", nothing proposed" : ""),
          color: colour ?? undefined,
          x: { domain, label: "Sampling iteration" },
          y: { domain: [0, 1], label: "acceptance", grid: true },
          marks: [
            Plot.frame({ stroke: "#d6d5d0" }),
            ...timeMarks(p, rules),
            Plot.rect([{ y1: band[0], y2: band[1] }], { x1: domain[0], x2: domain[1], y1: "y1", y2: "y2", fill: BAND, fillOpacity: 0.1 }),
            Plot.line(data, { x: "x", y: "rate", z: "target", stroke, strokeWidth: 1.5, strokeOpacity: colour ? 1 : 0.35 }),
            Plot.dot(data, { x: "x", y: "rate", fill: stroke, r: 2.5, fillOpacity: colour ? 1 : 0.5 }),
            Plot.tip(data, Plot.pointer({
              x: "x",
              y: "rate",
              title: (r: (typeof data)[number]) =>
                `${r.target}\n${range(r)}\naccepted ${r.accepted} of ${r.proposed}` + (r.rate === null ? "" : ` (${fmt(r.rate)})`),
            })),
          ],
        });
      }),
    );
  }

  function drawSwaps(p: SamplerPanel, rules: EventRule[]) {
    if (!p.swaps.length) {
      $(".swaps").textContent = "This Chain ran one temperature slot: no swaps.";
      return;
    }
    const data = p.swaps.map((s) => ({ ...s, x: windowX(s), pair: `slots ${s.temperature_pair}–${s.temperature_pair + 1}` }));
    const pairs = [...new Set(data.map((d) => d.pair))];
    const colour = colours(pairs);
    $(".swaps").replaceChildren(
      Plot.plot({
        width: 640,
        height: 220,
        marginLeft: 44,
        marginTop: 24,
        color: colour ?? undefined,
        x: { domain: timeDomain(p.windows), label: "Sampling iteration" },
        y: { domain: [0, 1], label: "swap acceptance", grid: true },
        marks: [
          Plot.frame({ stroke: "#d6d5d0" }),
          ...timeMarks(p, rules),
          Plot.line(data, { x: "x", y: "rate", z: "pair", stroke: colour ? "pair" : MUTED, strokeWidth: 1.5 }),
          Plot.dot(data, { x: "x", y: "rate", fill: colour ? "pair" : MUTED, r: 3 }),
          Plot.tip(data, Plot.pointer({
            x: "x",
            y: "rate",
            title: (s: (typeof data)[number]) =>
              `Temperature ${s.pair}\n${range(s)}\naccepted ${s.accepted} of ${s.attempted}` + (s.rate === null ? "" : ` (${fmt(s.rate)})`),
          })),
        ],
      }),
    );
  }

  function drawLadder(p: SamplerPanel, rules: EventRule[]) {
    const points = ladderPoints(p.ladder, p.burn_in_ended).map((d) => ({ ...d, name: `slot ${d.slot}` }));
    const skipped = p.ladder.filter((s) => !s.applied).length;
    const note = skipped ? ` ${skipped} adaptation${skipped === 1 ? "" : "s"} not applied (the new ladder had NaNs).` : "";
    if (!points.length) {
      $(".ladder").textContent = `This Chain ran one temperature slot: no ladder.${note}`;
      return;
    }
    const slots = [...new Set(points.map((d) => d.name))];
    const colour = colours(slots);
    const final = p.burn_in_ended?.temperatures ?? [];
    const domain = timeDomain(p.windows);
    // The final ladder holds for the rest of the Run: its line runs on to the end of the time axis.
    const line = [...points, ...points.filter((d) => d.final).map((d) => ({ ...d, iteration: domain[1] }))];
    $(".ladder").replaceChildren(
      Plot.plot({
        width: 640,
        height: 220,
        marginLeft: 44,
        marginTop: 24,
        color: colour ?? undefined,
        x: { domain, label: "Sampling iteration" },
        y: { domain: [0, 1], label: "temperature", grid: true },
        marks: [
          Plot.frame({ stroke: "#d6d5d0" }),
          ...timeMarks(p, rules),
          Plot.line(line, { x: "iteration", y: "temperature", z: "name", stroke: colour ? "name" : MUTED, strokeWidth: 1.5, curve: "step-after" }),
          Plot.dot(points, { x: "iteration", y: "temperature", fill: colour ? "name" : MUTED, r: 3 }),
          Plot.tip(points, Plot.pointer({
            x: "iteration",
            y: "temperature",
            title: (d: (typeof points)[number]) =>
              `Temperature ${d.name}: ${fmt(d.temperature)}\n` +
              (d.final ? `final ladder, burn-in end at iteration ${d.iteration}` : `adapted at iteration ${d.iteration}`) +
              (d.barrier === null ? "" : `\nswap barrier to slot ${d.slot + 1}: ${fmt(d.barrier)}`),
          })),
        ],
      }),
      el("p", {}, (final.length ? `Final ladder: ${final.map(fmt).join(", ")}.` : "") + note),
    );
  }

  function drawVariances(p: SamplerPanel) {
    const end = p.burn_in_ended;
    $(".variances-title").textContent = end
      ? `Final proposal variances (burn-in end, iteration ${end.iteration})`
      : "Final proposal variances";
    const rows = varianceRows(end);
    $(".variances").replaceChildren(
      !end
        ? "Burn-in has not ended in a closed Window."
        : !rows.length
          ? "No sampler has a proposal scale."
          : table(
              ["Sampler type", "Target", "Scales", "Proposal variances"],
              rows.map((r) => [r.type, r.target, String(r.variances.length), r.variances.map(fmt).join(", ")]),
              { numeric: [2] },
            ),
    );
  }

  function drawEvents(p: SamplerPanel) {
    const key = $(".events-key");
    const types = [...new Set(p.events.map((e) => e.type))];
    key.replaceChildren(
      el("span", { className: "muted" }, "Run events, a rule on every time axis (hover for details): "),
      ...types.map((type) => {
        const severity = severityOf(type);
        const span = el("span", { className: "event-key" });
        const swatch = el("i", { className: severity });
        swatch.style.borderColor = SEVERITY_COLOUR[severity];
        span.append(swatch, `${EVENT_NAMES[type] ?? type} (${p.events.filter((e) => e.type === type).length})`);
        return span;
      }),
    );
    $(".events").replaceChildren(
      p.events.length
        ? table(
            ["Iteration", "Phase", "Window", "Event", "Wall time (UTC)"],
            p.events.map((e) => [String(e.iteration), e.phase.replace("_", "-"), String(e.window), describeEvent(e), e.wall_time]),
            { numeric: [0, 2] },
          )
        : "No non-finite values, interrupts, burn-in end, resumes or Run end in the closed Windows.",
    );
  }

  function draw() {
    if (!panel) return;
    const rules = eventRules(panel.events);
    drawEvents(panel);
    drawAcceptance();
    drawSwaps(panel, rules);
    drawLadder(panel, rules);
    drawVariances(panel);
  }

  function load() {
    const mine = ++request;
    getJson<SamplerPanel>("/api/sampler", { chain })
      .then((p) => {
        if (mine !== request) return;
        panel = p;
        draw();
      })
      .catch((e: Error) => ($(".acceptance").textContent = `Could not load the sampler panel: ${e.message}`));
  }

  load();
};

export default render;
