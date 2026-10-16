// A pair of Infections A and B: P(A before B) per Chain with switch counts and the order per Draw; both Edge
// marginals per Chain; co-inclusion of A and B in third Infections' Parent sets; infection-time intervals and the
// lag; genotypes per Locus side by side with the differing Loci counted and highlighted. The Pair Inspector shows it
// by tab; the Pairs workspace shows all of it as its main content.

import * as Plot from "@observablehq/plot";
import { useCallback, useState, type ReactNode } from "react";
import { usePair } from "../../api/pairs";
import { useTheme } from "../../app/theme";
import { EntityLink } from "../../components/EntityLink";
import { PanelState } from "../../components/Panel";
import { PlotFigure, plotStyle, type PlotContext } from "../../components/PlotFigure";
import { formatFixed, formatProbability, formatSpan } from "../../lib/format";
import { differingLoci, type PairResponse } from "../../lib/pair";
import { chainColour, chainScale } from "../../lib/palette";
import { ChainLabel, Section } from "../layout";

const time = (x: number | null) => formatFixed(x, 1);

export type PairPart = "summary" | "genotypes" | "all";

function Table({ label, head, children }: { label: string; head: ReactNode[]; children: ReactNode }) {
  return (
    <table className="w-full text-table tabular-nums" aria-label={label}>
      <thead className="text-ink-2">
        <tr className="border-b border-line">
          {head.map((h, i) => (
            <th key={i} className={`py-0.5 font-medium ${i ? "text-right" : "text-left"}`}>{h}</th>
          ))}
        </tr>
      </thead>
      <tbody>{children}</tbody>
    </table>
  );
}

const td = (i: number) => (i ? "num" : "py-0.5");

function OrderSection({ p }: { p: PairResponse }) {
  const { a, b } = p;
  const { chart } = useTheme();
  const chains = p.order.chains.map((c) => c.chain);
  const draw = useCallback(
    ({ width, theme }: PlotContext) => {
      const points = p.order.chains.flatMap((c) => c.iteration.map((iteration, i) => ({ chain: c.chain, iteration, before: c.before[i] ?? 0 })));
      const label = (d: number) => (d ? `${a} first` : `${b} first`);
      return Plot.plot({
        width,
        height: 44 * chains.length + 28,
        // Room for "<id> first" on the left and each Chain's label on the right, at 11 px.
        marginLeft: Math.max(64, 6.5 * Math.max(a.length, b.length) + 40),
        marginRight: Math.max(16, 6.5 * Math.max(...chains.map((c) => c.length)) + 12),
        marginTop: 4,
        marginBottom: 26,
        style: plotStyle(theme),
        color: chainScale(theme, chains),
        x: { label: "Sampling iteration →", labelAnchor: "right", grid: true, tickFormat: "~s" },
        y: { domain: [0, 1], inset: 6, ticks: [0, 1], tickFormat: label, label: null },
        fy: { label: null, domain: chains, padding: 0.25, axis: "right", tickSize: 0, tickPadding: 6 },
        marks: [
          Plot.frame({ stroke: theme.grid }),
          Plot.lineY(points, { x: "iteration", y: "before", fy: "chain", stroke: "chain", strokeWidth: 2, curve: "step-after" }),
          Plot.tip(points, Plot.pointerX({ x: "iteration", y: "before", fy: "chain", channels: { Chain: "chain" }, format: { y: label, fy: false, stroke: false } })),
        ],
      });
    },
    [p, a, b, chains],
  );
  const draws = p.order.chains.reduce((n, c) => n + c.draws, 0);
  const switches = p.order.chains.reduce((n, c) => n + c.switches, 0);
  return (
    <Section title="Order" aside={`P(${a} before ${b})`}>
      <Table label={`P(${a} before ${b}) per Chain`} head={["Chain", "Draws", `P(${a} before ${b})`, "Switches"]}>
        {p.order.chains.map((c, k) => (
          <tr key={c.chain} className="border-b border-line/60">
            <td className={td(0)}><ChainLabel chain={c.chain} colour={chainColour(chart, k)} /></td>
            <td className={td(1)}>{c.draws}</td>
            <td className={td(2)}>{formatProbability(c.probability)}</td>
            <td className={td(3)}>{c.switches}</td>
          </tr>
        ))}
        <tr className="font-medium">
          <td className={td(0)}>All Chains</td>
          <td className={td(1)}>{draws}</td>
          <td className={td(2)}>{formatProbability(p.order.probability)}</td>
          <td className={td(3)}>{switches}</td>
        </tr>
      </Table>
      <p className="mt-1 text-axis text-muted">Switches: adjacent Draws at which the order changes. The line is up at Draws where {a} precedes {b} in the Ordering.</p>
      <PlotFigure plot={draw} height={44 * chains.length + 28} fallbackWidth={420} label={`Order of ${a} and ${b} per Draw, per Chain`} />
    </Section>
  );
}

function EdgesSection({ p }: { p: PairResponse }) {
  const { a, b } = p;
  const { chart } = useTheme();
  const ab = <EntityLink entity={{ kind: "edge", parent: a, child: b }} />;
  const ba = <EntityLink entity={{ kind: "edge", parent: b, child: a }} />;
  return (
    <Section title="Edge marginals">
      <Table label="Edge marginals per Chain" head={["Chain", ab, ba]}>
        {p.edges.chains.map((c, k) => (
          <tr key={c.chain} className="border-b border-line/60">
            <td className={td(0)}><ChainLabel chain={c.chain} colour={chainColour(chart, k)} /></td>
            <td className={td(1)}>{formatProbability(c.a_to_b)}</td>
            <td className={td(2)}>{formatProbability(c.b_to_a)}</td>
          </tr>
        ))}
        <tr className="font-medium">
          <td className={td(0)}>All Chains</td>
          <td className={td(1)}>{formatProbability(p.edges.a_to_b)}</td>
          <td className={td(2)}>{formatProbability(p.edges.b_to_a)}</td>
        </tr>
      </Table>
      <p className="mt-1 text-axis text-muted">Mean over Draws of the summed weight of the child's Parent sets holding the parent.</p>
    </Section>
  );
}

function CoInclusionSection({ p }: { p: PairResponse }) {
  const { a, b } = p;
  const chains = p.order.chains.map((c) => c.chain);
  return (
    <Section title="Co-inclusion" aside={`of ${a} and ${b} in a third Infection's Parent sets`}>
      {p.co_inclusion.length ? (
        <Table label="Co-inclusion" head={["Infection", "All Chains", ...chains]}>
          {p.co_inclusion.map((r) => (
            <tr key={r.infection} className="border-b border-line/60">
              <td className={td(0)}><EntityLink entity={{ kind: "infection", id: r.infection }} /></td>
              <td className={td(1)}>{formatProbability(r.weight)}</td>
              {chains.map((ch, k) => (
                <td key={ch} className={td(k + 2)}>{formatProbability(r.chains.find((c) => c.chain === ch)?.weight ?? 0)}</td>
              ))}
            </tr>
          ))}
        </Table>
      ) : (
        <p className="text-ink-2">No Infection has a Parent set holding both {a} and {b}.</p>
      )}
    </Section>
  );
}

const pairTimes = (p: PairResponse) => [{ infection: p.a, ...p.times.a }, { infection: p.b, ...p.times.b }];

function TimesSection({ p }: { p: PairResponse }) {
  const { a, b } = p;
  const times = pairTimes(p);
  const draw = useCallback(
    ({ width, theme }: PlotContext) => {
      const times = pairTimes(p);
      const intervals = times.filter((t) => t.lower !== null);
      return Plot.plot({
        width,
        height: 92,
        marginLeft: 64,
        marginRight: 12,
        marginBottom: 28,
        style: plotStyle(theme),
        x: { label: "Time →", labelAnchor: "right", grid: true },
        y: { domain: [a, b], label: null },
        marks: [
          Plot.link(intervals, { x1: "lower", x2: "upper", y1: "infection", y2: "infection", stroke: theme.ink2, strokeWidth: 2, strokeLinecap: "round" }),
          Plot.dot(intervals, { x: "median", y: "infection", fill: theme.ink2, r: 4, stroke: theme.surface, strokeWidth: 1.5 }),
          Plot.dot(times, { x: "observation_time", y: "infection", symbol: "times", stroke: theme.ink, r: 5, strokeWidth: 2 }),
          Plot.tip(times, Plot.pointerY({ x: "median", y: "infection", channels: { observed: "observation_time", lower: "lower", upper: "upper" }, format: { x: ".1f", y: true } })),
        ],
      });
    },
    [p, a, b],
  );
  return (
    <Section title="Infection times">
      <Table label="Infection times" head={["Infection", "Observed", "Infection time, median", "95% interval"]}>
        {times.map((t) => (
          <tr key={t.infection} className="border-b border-line/60">
            <td className={td(0)}><EntityLink entity={{ kind: "infection", id: t.infection }} /></td>
            <td className={td(1)}>{time(t.observation_time)}</td>
            <td className={td(2)}>{time(t.median)}</td>
            <td className={td(3)}>{formatSpan(t.lower, t.upper, time)}</td>
          </tr>
        ))}
      </Table>
      <p className="mt-1">Lag (observation time of {b} minus {a}): <strong className="tabular-nums">{time(p.times.lag)}</strong></p>
      <PlotFigure plot={draw} height={92} fallbackWidth={420} label={`Infection-time intervals of ${a} and ${b}`} />
      <p className="text-axis text-muted">Line: 95% posterior interval of infection time; dot: median; ×: observation time.</p>
    </Section>
  );
}

/** The sticky Locus cell of a differing row: the row's warning tint over the opaque surface, and its left bar. */
const DIFFERS_STICKY = {
  background: "linear-gradient(color-mix(in srgb, var(--warning) 15%, transparent), color-mix(in srgb, var(--warning) 15%, transparent)), var(--surface)",
  boxShadow: "inset 2px 0 0 var(--warning), inset -1px 0 0 var(--line)",
};

/** A genotype's Allele bits; `missing` where there is none ("unobserved" only for an Observed genotype). */
function Bits({ bits, missing = "–" }: { bits: string | null | undefined; missing?: string }) {
  return bits ? <span className="id">{bits}</span> : <span className="text-muted">{missing}</span>;
}

function GenotypesSection({ p }: { p: PairResponse }) {
  const { a, b } = p;
  const [onlyDiffering, setOnlyDiffering] = useState(false);
  const differing = differingLoci(p.genotypes.rows);
  const truth = p.genotypes.truth;
  const shown = onlyDiffering ? p.genotypes.rows.filter((r) => r.differs) : p.genotypes.rows;
  const head = ["Locus", `Observed ${a}`, `Latent ${a}`, "Agr.", `Latent ${b}`, "Agr.", `Observed ${b}`, ...(truth ? [`True ${a}`, `True ${b}`] : [])];
  return (
    <Section
      title="Genotypes"
      aside={
        <label className="flex items-center gap-1 text-ink-2">
          <input type="checkbox" checked={onlyDiffering} onChange={(e) => setOnlyDiffering(e.target.checked)} className="accent-[var(--accent)]" />
          Differing only
        </label>
      }
    >
      <p className={`mb-1 ${differing.count ? "font-medium text-warning-ink" : "text-ink-2"}`}>{differing.text}.</p>
      <p className="mb-1 text-axis text-muted">Latent: posterior mode over every Chain's Draws; Agr.: the fraction of Draws that had it.</p>
      {/* Wider than the Inspector: it scrolls sideways here, the Locus column staying in view (as DataTable's). */}
      <div className="-mx-3 overflow-x-auto">
        <table className="min-w-full text-table tabular-nums" aria-label="Genotypes per Locus">
          <thead className="text-ink-2">
            <tr className="border-b border-line">
              {head.map((h, i) => (
                <th key={i} className={`whitespace-nowrap py-0.5 pr-2 text-left font-medium ${h === "Agr." ? "text-right" : ""} ${i === 0 ? "dt-sticky pl-3" : ""}`}>{h}</th>
              ))}
            </tr>
          </thead>
          <tbody>
            {shown.map((r) => (
              <tr
                key={r.locus}
                data-differs={r.differs}
                title={r.differs ? "The Latent modes differ" : undefined}
                className={`border-b border-line/60 ${r.differs ? "bg-warning/15 font-medium" : ""}`}
              >
                <td className="dt-sticky whitespace-nowrap py-0.5 pl-3 pr-2" style={r.differs ? DIFFERS_STICKY : undefined}>{r.locus}</td>
                <td className="pr-2"><Bits bits={r.observed_a} missing="unobserved" /></td>
                <td className="pr-2"><Bits bits={r.latent_a} /></td>
                <td className="num pr-2">{formatProbability(r.agreement_a)}</td>
                <td className="pr-2"><Bits bits={r.latent_b} /></td>
                <td className="num pr-2">{formatProbability(r.agreement_b)}</td>
                <td className="pr-2"><Bits bits={r.observed_b} missing="unobserved" /></td>
                {truth && <td className="pr-2"><Bits bits={r.true_a} /></td>}
                {truth && <td className="pr-3"><Bits bits={r.true_b} /></td>}
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </Section>
  );
}

export function PairView({ a, b, part = "all" }: { a: string; b: string; part?: PairPart }) {
  const { data: p, error, isLoading } = usePair(a, b);
  if (error || isLoading || !p) return <PanelState loading={isLoading} error={error} empty="Pick two different Infections." />;
  return (
    <div className="text-ui">
      {part !== "genotypes" && (
        <>
          <OrderSection p={p} />
          <EdgesSection p={p} />
          <CoInclusionSection p={p} />
          <TimesSection p={p} />
        </>
      )}
      {part !== "summary" && <GenotypesSection p={p} />}
    </div>
  );
}
