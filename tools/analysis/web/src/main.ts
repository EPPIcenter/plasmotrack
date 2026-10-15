// The page layout: a header with the Analysis summary and the navigation bar, and the current view below it.

import "./style.css";
import { getAnalysis, type AnalysisSummary } from "./api";
import { viewFor, views } from "./registry";
import { formatRoute, parseRoute } from "./router";

const header = document.querySelector<HTMLElement>("#header")!;
const main = document.querySelector<HTMLElement>("#view")!;
let cleanup: void | (() => void);

function drawHeader(analysis: AnalysisSummary, current: string) {
  const counts = `${analysis.chains.length} Chains, ${analysis.draws.toLocaleString()} Draws, ${analysis.infections.toLocaleString()} Infections, ${analysis.loci.toLocaleString()} Loci`;
  header.innerHTML = `<h1></h1><span class="muted"></span><nav></nav>`;
  header.querySelector("h1")!.textContent = analysis.label;
  header.querySelector("span")!.textContent = counts;
  header.querySelector("nav")!.append(
    ...views
      .filter((v) => v.nav)
      .map((v) => Object.assign(document.createElement("a"), {
        href: formatRoute(v.route),
        textContent: v.title,
        className: v.route === current ? "current" : "",
      })),
  );
}

async function show(analysis: AnalysisSummary) {
  const { route, params } = parseRoute(location.hash);
  const entry = viewFor(route);
  drawHeader(analysis, entry.route);
  document.title = `${entry.title} · ${analysis.label}`;
  if (typeof cleanup === "function") cleanup();
  const { default: render } = await entry.load();
  main.replaceChildren();
  cleanup = render(main, {
    params,
    analysis,
    navigate: (to, p) => (location.hash = formatRoute(to, p)),
    replaceParams: (p) => history.replaceState(null, "", formatRoute(entry.route, p)),
  });
}

getAnalysis()
  .then((analysis) => {
    addEventListener("hashchange", () => show(analysis));
    return show(analysis);
  })
  .catch((e: Error) => (main.textContent = `Could not reach the pt-analysis backend: ${e.message}`));
