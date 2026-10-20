// Genotypes of the selected Infections (the network's multi-selection, `ids`): one column per Infection, one row per
// Locus, each cell an allele strip. Opens on the Loci where the Latent modes are not all the same.
//
// URL: ?ids=a,b,c&loci=all

import { useState, type FormEvent } from "react";
import { useGroupGenotypes } from "../../api/genotypes";
import { useNetwork } from "../../api/network";
import { useSearchState } from "../../app/search";
import { useMultiSelection } from "../../app/selection";
import { AlleleStrip } from "../../components/AlleleStrip";
import { EntityLink } from "../../components/EntityLink";
import { Panel, PanelState } from "../../components/Panel";
import { WorkspaceLink } from "../../components/WorkspaceLink";
import { formatCount } from "../../lib/format";

const MAX = 12;

export default function GenotypesWorkspace() {
  const selection = useMultiSelection();
  const [search, set] = useSearchState();
  const ids = selection.ids.slice(0, MAX);
  const genotypes = useGroupGenotypes(ids);
  const network = useNetwork();
  const [text, setText] = useState("");
  const onlyDiffering = search.loci !== "all";
  const loci = genotypes.data?.loci ?? [];
  const shown = onlyDiffering ? loci.filter((l) => l.differs) : loci;

  const add = (e: FormEvent) => {
    e.preventDefault();
    const id = text.trim();
    if (!id || selection.ids.includes(id) || selection.ids.length >= MAX) return;
    selection.set([...selection.ids, id]);
    setText("");
  };

  return (
    <div className="h-full p-2">
      <Panel
        className="h-full"
        title="Genotypes"
        description="Filled slot: an Allele of the Latent mode. Outline: the Observed genotype has it. The same slot is the same Allele."
        noExport
        bodyClassName="flex flex-col"
        actions={
          <label className="flex items-center gap-1 text-axis text-ink-2">
            <input type="checkbox" checked={onlyDiffering} onChange={(e) => set({ loci: e.target.checked ? undefined : "all" })} className="accent-[var(--accent)]" />
            Differing only
          </label>
        }
      >
        <div className="flex shrink-0 flex-wrap items-center gap-2 border-b border-line px-2 py-1 text-axis">
          <form onSubmit={add} className="flex items-center gap-1">
            <input
              aria-label="Add an Infection"
              list="genotype-ids"
              value={text}
              onChange={(e) => setText(e.target.value)}
              placeholder="Add Infection"
              className="h-5 w-36 rounded-[3px] border border-line bg-surface px-1.5"
            />
            <datalist id="genotype-ids">
              {(network.data?.nodes ?? []).map((n) => <option key={n.id} value={n.id} />)}
            </datalist>
          </form>
          <WorkspaceLink to="network" className="text-accent hover:underline">Choose on the Network</WorkspaceLink>
          {selection.ids.length > MAX && <span className="text-warning-ink">The first {MAX} Infections are shown.</span>}
        </div>
        {ids.length === 0 ? (
          <p className="px-2 py-3 text-axis text-muted">Select Infections on the Network, or name one above. Shift-click nodes to add them.</p>
        ) : genotypes.error || genotypes.isLoading || !genotypes.data ? (
          <PanelState loading={genotypes.isLoading} error={genotypes.error} />
        ) : (
          <>
            <p className={`shrink-0 px-2 py-1 text-axis ${genotypes.data.differing ? "font-medium text-warning-ink" : "text-ink-2"}`}>
              {genotypes.data.differing
                ? `${formatCount(genotypes.data.differing)} of ${formatCount(loci.length)} Loci differ across ${formatCount(ids.length)} Infections.`
                : `No Locus of ${formatCount(loci.length)} differs across ${formatCount(ids.length)} Infections.`}
            </p>
            {shown.length === 0 ? (
              <p className="px-2 py-2 text-axis text-muted">No Locus where the Latent modes differ.</p>
            ) : (
              <div className="min-h-0 flex-1 overflow-auto">
                <table className="text-table" aria-label="Genotypes of the selection">
                  <thead className="sticky top-0 z-10 bg-surface text-ink-2">
                    <tr className="border-b border-line">
                      <th className="dt-sticky bg-surface py-1 pl-2 pr-3 text-left font-medium">Locus</th>
                      {ids.map((id) => (
                        <th key={id} className="whitespace-nowrap py-1 pr-3 text-left font-medium">
                          <EntityLink entity={{ kind: "infection", id }} />
                          <button type="button" aria-label={`Remove ${id}`} onClick={() => selection.toggle(id)} className="ml-1 text-muted hover:text-ink">×</button>
                        </th>
                      ))}
                    </tr>
                  </thead>
                  <tbody>
                    {shown.map((locus) => (
                      <tr key={locus.locus} data-differs={locus.differs} className={`border-b border-line/60 ${locus.differs ? "bg-warning/15" : ""}`}>
                        <td className="dt-sticky max-w-52 truncate whitespace-nowrap bg-surface py-0.5 pl-2 pr-3" title={locus.locus}>{locus.locus}</td>
                        {locus.members.map((member) => (
                          <td key={member.id} className="py-0.5 pr-3">
                            <AlleleStrip latent={member.mode} observed={member.observed} agreement={member.agreement} alleles={locus.alleles} missing="unobserved" />
                          </td>
                        ))}
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            )}
          </>
        )}
      </Panel>
    </div>
  );
}
