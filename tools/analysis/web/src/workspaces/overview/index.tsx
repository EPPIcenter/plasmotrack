// Overview: the health check and the headline posterior on one screen. KPIs; the log posterior's trace and rank
// histograms; ESS and R-hat strips per convergence group, each opening Convergence filtered to its group; the
// attention list (worst-mixing quantities, edges and Parent sets the Chains disagree on, sampler types outside the
// band, Run events); and the headline posterior (composition, q and effective parent-set size over Infections, top
// edges, a network thumbnail).

import { useOverview } from "../../api/overview";
import { Panel, PanelState } from "../../components/Panel";
import { QuantityTrace } from "../../components/QuantityTrace";
import { overviewTables } from "../../lib/overview";
import { Attention } from "./Attention";
import { Headline } from "./Headline";
import { Kpis } from "./Kpis";
import { MixingByGroup } from "./MixingByGroup";

const LOG_POSTERIOR = { kind: "quantity" as const, group: "log posterior", name: "log posterior" };

export default function OverviewWorkspace() {
  const { data: o, error, isLoading } = useOverview();
  return (
    <div className="@container flex min-h-full min-w-0 flex-col gap-2 p-2">
      <Kpis o={o} />
      {error || isLoading || !o ? (
        <PanelState loading={isLoading} error={error} />
      ) : (
        <>
          <div className="grid grid-cols-1 gap-2 @5xl:grid-cols-[minmax(0,5fr)_minmax(0,6fr)]">
            <Panel title="Log posterior" description="Trace per Chain over Sampling iterations; brush for ESS and R-hat of a range">
              <QuantityTrace entity={LOG_POSTERIOR} height={150} compact />
            </Panel>
            <Panel title="Mixing by group" description="Every quantity's pooled ESS and R-hat; accent tick: median ESS; rules: amber and red thresholds; a row opens Convergence" data={() => overviewTables(o)} exportName="mixing">
              <MixingByGroup o={o} />
            </Panel>
          </div>
          <Panel title="Attention" description="Rows select into the Inspector" data={() => overviewTables(o)} exportName="attention">
            <Attention o={o} />
          </Panel>
          <Panel title="Posterior headline" description="Every Infection pooled over Chains and Draws" data={() => overviewTables(o)} exportName="headline">
            <Headline o={o} />
          </Panel>
        </>
      )}
    </div>
  );
}
