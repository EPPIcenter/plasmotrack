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
