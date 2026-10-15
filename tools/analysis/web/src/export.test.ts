import { describe, expect, it } from "vitest";
import { composeSvg, exportFileName, htmlTable, intervalColumns, intervalFields, recordsTable, toCsv } from "./export";

describe("toCsv", () => {
  it("writes a header and one line per row, quoting only what needs it", () => {
    const csv = toCsv({ columns: ["child", "parent", "note"], rows: [["9", "1", 'says "hi", twice'], ["9", "2", "two\nlines"]] });
    expect(csv).toBe('child,parent,note\r\n9,1,"says ""hi"", twice"\r\n9,2,"two\nlines"\r\n');
  });

  it("writes missing values empty, non-finite numbers as R reads them, booleans as true and false", () => {
    const csv = toCsv({ columns: ["a", "b", "c", "d", "e", "f"], rows: [[null, undefined, Infinity, -Infinity, NaN, true], [0.125, 3, "", false, "x", 1e-7]] });
    expect(csv).toBe("a,b,c,d,e,f\r\n,,Inf,-Inf,NaN,true\r\n0.125,3,,false,x,1e-7\r\n");
  });
});

describe("recordsTable", () => {
  it("lays records out in the given column order, a missing key as an empty cell", () => {
    const t = recordsTable(["group", "name", "chain", "ess"], [{ name: "q", group: "g", ess: 4, chain: null }, { group: "g", name: "r" }]);
    expect(toCsv(t)).toBe("group,name,chain,ess\r\ng,q,,4\r\ng,r,,\r\n");
  });
});

describe("interval columns", () => {
  it("spread an interval over three columns named after its field, empty when there is none", () => {
    expect(intervalColumns("psi_edge")).toEqual(["psi_edge_mean", "psi_edge_lower", "psi_edge_upper"]);
    expect(intervalFields("psi_edge", { mean: 0.5, lower: 0.25, upper: 0.75 })).toEqual({ psi_edge_mean: 0.5, psi_edge_lower: 0.25, psi_edge_upper: 0.75 });
    expect(intervalFields("mean_q", null)).toEqual({ mean_q_mean: null, mean_q_lower: null, mean_q_upper: null });
  });
});

describe("exportFileName", () => {
  it("joins the Analysis label, the route and the query state", () => {
    const params = new URLSearchParams("child=9&sort=mean:chain A&desc=1");
    expect(exportFileName("Golden fixture", "edges", params, "csv")).toBe("Golden-fixture_edges_child=9_sort=mean-chain-A_desc=1.csv");
  });

  it("leaves out an empty query and appends a part name when a view gives several files", () => {
    expect(exportFileName("golden", "convergence", new URLSearchParams(), "svg")).toBe("golden_convergence.svg");
    expect(exportFileName("golden", "pair", new URLSearchParams("a=1&b=2"), "csv", "Genotypes per Locus"))
      .toBe("golden_pair_a=1_b=2_Genotypes-per-Locus.csv");
  });

  it("keeps characters no file system allows out of the name, and the name short", () => {
    const name = exportFileName('a/b\\c:d*e?"f<g>h|i', "network", new URLSearchParams(`q=${"x".repeat(400)}`), "png");
    expect(name).toMatch(/^a-b-c-d-e-f-g-h-i_network_q=x+\.png$/);
    expect(name.length).toBeLessThanOrEqual(200);
  });
});

// The shape of an HTMLTableElement that htmlTable reads: its rows' cells, their tag and text.
const cell = (tagName: "TH" | "TD", textContent: string) => ({ tagName, textContent });
const tr = (...cells: ReturnType<typeof cell>[]) => ({ cells });

describe("htmlTable", () => {
  it("takes a header row of th cells as the column names and each later row's text as it was shown", () => {
    const t = htmlTable({ rows: [tr(cell("TH", "Chain"), cell("TH", " q ")), tr(cell("TD", "A"), cell("TD", "0.512\n  ")), tr(cell("TD", "B"), cell("TD", "–"))] });
    expect(toCsv(t)).toBe("Chain,q\r\nA,0.512\r\nB,–\r\n");
  });

  it("numbers the columns of a table without a header row and pads short rows", () => {
    const t = htmlTable({ rows: [tr(cell("TD", "a"), cell("TD", "b")), tr(cell("TD", "c"))] });
    expect(t.columns).toEqual(["1", "2"]);
    expect(t.rows).toEqual([["a", "b"], ["c", ""]]);
  });
});

describe("composeSvg", () => {
  const plot = '<svg xmlns="http://www.w3.org/2000/svg" width="300" height="100"><circle r="2"/></svg>';

  it("places each piece where it was on screen, cropped to the box around them all", () => {
    const composed = composeSvg([
      { kind: "svg", markup: plot, box: { x: 20, y: 50, width: 300, height: 100 } },
      { kind: "svg", markup: plot.replace("circle", "rect"), box: { x: 340, y: 50.4, width: 300, height: 100 } },
    ])!;
    expect([composed.width, composed.height]).toEqual([620, 100.4]);
    const svg = composed.markup;
    expect(svg).toMatch(/^<svg xmlns="http:\/\/www.w3.org\/2000\/svg" width="620" height="100.4" viewBox="0 0 620 100.4">/);
    expect(svg).toContain(`<g transform="translate(0,0)">${plot}</g>`);
    expect(svg).toContain('<g transform="translate(320,0.4)"><svg xmlns="http://www.w3.org/2000/svg" width="300" height="100"><rect');
    expect(svg.endsWith("</svg>")).toBe(true);
  });

  it("draws a canvas as an image and HTML as a foreign object, in the order given", () => {
    const svg = composeSvg([
      { kind: "image", href: 'data:image/png;base64,AA&"', box: { x: 0, y: 0, width: 10, height: 5 } },
      { kind: "html", markup: '<div xmlns="http://www.w3.org/1999/xhtml">q</div>', box: { x: 0, y: 5, width: 10, height: 5 } },
    ])!.markup;
    expect(svg).toContain('<image x="0" y="0" width="10" height="5" href="data:image/png;base64,AA&amp;&quot;"/>');
    expect(svg).toContain('<foreignObject x="0" y="5" width="10" height="5"><div xmlns="http://www.w3.org/1999/xhtml">q</div></foreignObject>');
    expect(svg.indexOf("<image")).toBeLessThan(svg.indexOf("<foreignObject"));
  });

  it("gives nothing when there is nothing to draw", () => {
    expect(composeSvg([])).toBeNull();
  });
});
