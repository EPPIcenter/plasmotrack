// The Inspector's frame: the entity's kind and label, "Open as page" and close, its tabs (from the entity registry),
// and the kind's module rendering the current tab. Selecting an entity never changes the workspace.

import { Tabs } from "radix-ui";
import { Component, lazy, Suspense, type ComponentType, type ReactNode } from "react";
import { useSelection } from "../app/selection";
import { useLinkTo } from "../app/search";
import { entityLabel, formatEntity, type Entity, type EntityKind } from "../lib/entity";
import { entities, type InspectorProps } from "./registry";

const loaded = new Map<EntityKind, ComponentType<InspectorProps>>();
function view(kind: EntityKind): ComponentType<InspectorProps> {
  if (!loaded.has(kind)) loaded.set(kind, lazy(entities[kind].load as () => Promise<{ default: ComponentType<InspectorProps> }>));
  return loaded.get(kind)!;
}

class Boundary extends Component<{ children: ReactNode; reset: string }, { error: Error | null; reset: string }> {
  state = { error: null as Error | null, reset: this.props.reset };
  static getDerivedStateFromError(error: Error) {
    return { error };
  }
  static getDerivedStateFromProps(props: { reset: string }, state: { reset: string }) {
    return props.reset !== state.reset ? { error: null, reset: props.reset } : null;
  }
  render() {
    return this.state.error ? <p className="p-3 text-axis text-critical-ink">Could not show this: {this.state.error.message}</p> : this.props.children;
  }
}

export function Inspector({ entity, page = false }: { entity: Entity; page?: boolean }) {
  const { tab, setTab, select } = useSelection();
  const { go } = useLinkTo();
  const entry = entities[entity.kind];
  const current = entry.tabs.find((t) => t.id === tab)?.id ?? entry.tabs[0]!.id;
  const View = view(entity.kind);
  const key = formatEntity(entity);

  const openAsPage = () => go("inspect");

  return (
    <aside className="flex h-full min-h-0 flex-col bg-surface" aria-label={`Inspector: ${entry.title} ${entityLabel(entity)}`}>
      <header className="flex h-8 shrink-0 items-center gap-2 border-b border-line px-3">
        <span className="text-[10px] font-semibold uppercase tracking-wider text-muted">{entry.title}</span>
        <h2 className={`min-w-0 truncate font-semibold ${entity.kind === "quantity" || entity.kind === "stratum" ? "" : "id"}`} title={entityLabel(entity)}>
          {entityLabel(entity)}
        </h2>
        <div className="ml-auto flex items-center gap-1">
          {!page && (
            <button type="button" onClick={openAsPage} className="h-6 rounded-[3px] px-1.5 text-axis text-ink-2 hover:bg-hover" title="Open as the workspace, full width">
              Open as page ⤢
            </button>
          )}
          {!page && (
            <button type="button" onClick={() => select(null)} className="h-6 w-6 rounded-[3px] text-ink-2 hover:bg-hover" aria-label="Close the Inspector" title="Close">
              ×
            </button>
          )}
        </div>
      </header>
      <Tabs.Root value={current} onValueChange={setTab} className="flex min-h-0 flex-1 flex-col">
        {entry.tabs.length > 1 && (
          <Tabs.List className="flex h-7 shrink-0 gap-3 border-b border-line px-3" aria-label="Inspector tabs">
            {entry.tabs.map((t) => (
              <Tabs.Trigger
                key={t.id}
                value={t.id}
                onClick={() => t.id !== current && setTab(t.id)}
                className="-mb-px border-b-2 border-transparent text-ui text-ink-2 hover:text-ink data-[state=active]:border-accent data-[state=active]:text-ink"
              >
                {t.title}
              </Tabs.Trigger>
            ))}
          </Tabs.List>
        )}
        <div className="min-h-0 flex-1 overflow-y-auto">
          <Boundary reset={`${key}/${current}`}>
            <Suspense fallback={<p className="p-3 text-muted">Loading…</p>}>
              <View entity={entity} tab={current} />
            </Suspense>
          </Boundary>
        </div>
      </Tabs.Root>
    </aside>
  );
}
