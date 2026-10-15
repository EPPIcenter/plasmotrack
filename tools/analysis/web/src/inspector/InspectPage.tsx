// "Open as page": the selected entity as the workspace, full width (for wide content such as a genotype heatmap).

import { useSelection } from "../app/selection";
import { Inspector } from "./Inspector";

export default function InspectPage() {
  const { entity } = useSelection();
  if (!entity) return <p className="p-6 text-muted">Nothing is selected. Search with ⌘K, or select an entity in any workspace.</p>;
  return <Inspector entity={entity} page />;
}
