// The page: top bar, left rail of workspaces, the workspace and the Inspector side by side (resizable), status bar.
// The Inspector opens when the URL selects an entity and never changes the workspace.

import { Outlet, useLocation } from "@tanstack/react-router";
import { Suspense, useEffect, useState } from "react";
import { Group, Panel, Separator } from "react-resizable-panels";
import { CommandPalette } from "../components/CommandPalette";
import { Inspector } from "../inspector/Inspector";
import { Rail } from "./Rail";
import { useSelection } from "./selection";
import { StatusBar } from "./StatusBar";
import { TopBar } from "./TopBar";

export function Shell() {
  const { entity } = useSelection();
  const { pathname } = useLocation();
  const [palette, setPalette] = useState(false);
  const page = pathname === "/inspect";

  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      if ((e.metaKey || e.ctrlKey) && e.key.toLowerCase() === "k") {
        e.preventDefault();
        setPalette((open) => !open);
      }
    };
    addEventListener("keydown", onKey);
    return () => removeEventListener("keydown", onKey);
  }, []);

  return (
    <div className="grid h-full grid-rows-[36px_minmax(0,1fr)_22px] bg-plane text-ink">
      <TopBar onSearch={() => setPalette(true)} />
      <div className="flex min-h-0">
        <Rail />
        <Group orientation="horizontal" className="min-w-0 flex-1" id="shell">
          <Panel id="workspace" minSize={360} className="min-w-0">
            <main className="h-full min-h-0 overflow-auto [contain:inline-size]" aria-label="Workspace">
              <Suspense fallback={<p className="p-3 text-muted">Loading…</p>}>
                <Outlet />
              </Suspense>
            </main>
          </Panel>
          {entity && !page && (
            <>
              <Separator className="w-px bg-line outline-none transition-colors hover:bg-accent data-[separator=active]:bg-accent" />
              <Panel id="inspector" defaultSize="34" minSize={300} maxSize="70" className="min-w-0">
                <Inspector entity={entity} />
              </Panel>
            </>
          )}
        </Group>
      </div>
      <StatusBar />
      <CommandPalette open={palette} onOpenChange={setPalette} />
    </div>
  );
}
