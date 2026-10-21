// "Open as page": the selected entity as the workspace, full width. Every tab shows at once in a responsive grid
// (an Infection's summary, parents and genotype heatmap side by side), and a single-tab entity lays its sections out
// in columns. A link goes back to the workspace it was opened from (`from=`, the Overview by default).

import { Suspense } from "react";
import { useSearchState } from "../app/search";
import { useSelection } from "../app/selection";
import { workspaceById } from "../app/workspaces";
import { PanelState } from "../components/Panel";
import { WorkspaceLink } from "../components/WorkspaceLink";
import { entityLabel, formatEntity } from "../lib/entity";
import { Boundary, View } from "./Inspector";
import { parseFrom } from "./page";
import { entities } from "./registry";

export default function InspectPage() {
  const { entity } = useSelection();
  const [search] = useSearchState();
  const from = parseFrom(search.from);
  const back = (
    <WorkspaceLink to={from} className="text-axis text-accent hover:underline">
      ← Back to {workspaceById(from)?.title ?? "Overview"}
    </WorkspaceLink>
  );
  if (!entity) {
    return (
      <div className="p-6">
        <p className="text-muted">Nothing is selected. Search with ⌘K, or select an entity in any workspace.</p>
        <p className="mt-2">{back}</p>
      </div>
    );
  }
  const entry = entities[entity.kind];
  const key = formatEntity(entity);
  const label = entityLabel(entity);
  return (
    <article className="flex h-full min-h-0 flex-col bg-surface" aria-label={`${entry.title} ${label}`}>
      <header className="flex h-8 shrink-0 items-center gap-3 border-b border-line px-3">
        {back}
        <span className="text-[10px] font-semibold uppercase tracking-wider text-muted">{entry.title}</span>
        <h2 className={`min-w-0 truncate font-semibold ${entity.kind === "quantity" || entity.kind === "stratum" ? "" : "id"}`} title={label}>
          {label}
        </h2>
      </header>
      <div className="min-h-0 flex-1 overflow-auto">
        {entry.tabs.length > 1 ? (
          <div className="grid grid-cols-[repeat(auto-fill,minmax(440px,1fr))] items-start">
            {entry.tabs.map((t) => (
              <section key={t.id} aria-label={t.title} className="min-w-0 border-r border-b border-line">
                <h2 className="border-b border-line bg-plane px-3 py-1 text-ui font-semibold">{t.title}</h2>
                <Boundary reset={`${key}/${t.id}`}>
                  <Suspense fallback={<PanelState loading />}>
                    <View entity={entity} tab={t.id} page />
                  </Suspense>
                </Boundary>
              </section>
            ))}
          </div>
        ) : (
          <Boundary reset={key}>
            <Suspense fallback={<PanelState loading />}>
              <View entity={entity} tab={entry.tabs[0]?.id ?? ""} page />
            </Suspense>
          </Boundary>
        )}
      </div>
    </article>
  );
}
