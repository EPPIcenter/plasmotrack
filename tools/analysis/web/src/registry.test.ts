import { describe, expect, it } from "vitest";
import { viewFor, views } from "./registry";

describe("the view registry", () => {
  it("gives every view a distinct route and a title", () => {
    const routes = views.map((v) => v.route);
    expect(new Set(routes).size).toBe(routes.length);
    for (const v of views) expect(v.title).not.toBe("");
  });

  it("has the convergence table first, as the default view, and the trace panel", () => {
    expect(viewFor("").route).toBe("convergence");
    expect(viewFor("no such view").route).toBe("convergence");
    expect(viewFor("trace").title).toBe("Trace");
  });

  it("loads each view as a module whose default export renders it", async () => {
    for (const v of views) expect(typeof (await v.load()).default).toBe("function");
  });
});

// Exports (U11) come from the layout, not the views: SVG and PNG draw any view (its SVG elements and canvases, else
// its HTML), and CSV reads a view's `<table>` elements or its ViewHooks.table. This checks each registered view's
// source gives the default mechanism something, and that a virtualised table, which keeps only the rows in view in
// the DOM, hands over its data through the hook.

describe("every registered view exports", () => {
  // Read at run time, not through import.meta.glob, whose expansion vitest caches across runs.
  const source = async (route: string): Promise<string> =>
    (await import(`./views/${route}.ts?raw`).catch(() => ({ default: "" }))).default;
  const hasHook = (src: string) => /\btable: \(\) =>/.test(src);
  const drawsGraphics = (src: string) => /\bPlot\.\w+\(|new Sigma\b|createElementNS\([^)]*"svg"|createElement\("canvas"\)/.test(src);
  const rendersTables = (src: string) =>
    /createElement\("table"\)|\bel\(\s*"table"|<table\b/.test(src) || (/import \{[^}]*\btable\b[^}]*\} from "\.\.\/dom"/.test(src) && /\btable\(/.test(src));

  it.each(views.map((v) => v.route))("%s, as src/views/<route>.ts", async (route) => {
    const src = await source(route);
    expect(src, `src/views/${route}.ts`).not.toBe("");
    expect(hasHook(src) || drawsGraphics(src) || rendersTables(src)).toBe(true);
    if (/\bvirtualList\s*(<[^>]*>)?\(/.test(src)) expect(hasHook(src), "a virtualised table needs ViewHooks.table").toBe(true);
  });
});
